
import os
import pandas as pd
import argparse
import pysam
from collections import Counter
from collections import defaultdict
import re

# for the isoquant conda environment


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process peptide and SO ID mapping files.")

    parser.add_argument(
        "--classification_txt",
        required=True,
        help="Path to the peptides file."
    )

    parser.add_argument(
        "--custom_gtf",
        required=True,
        help="Path to the SO ID mapping file."
    )

    parser.add_argument(
        "--bam_file",
        required=True,
        help="Path to the SO ID mapping file."
    )

    parser.add_argument(
        "--isoquant_transcript_counts",
        required=True,
        help="Cell type name."
    )

    parser.add_argument(
        "--cell_type",
        required=True,
        help="Cell type name."
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Outdirectory for results."
    )

    return parser.parse_args()


################################################################################
# PATH DEFINITIONS
################################################################################
classification_txt = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/SQANTI3_QC/HUVEC/isoforms_classification.txt'
cell_type = 'HUVEC'
bam_file = '/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment/HUVEC/minimap2_align/merged/HUVEC_merged_sorted.bam'
outdir = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/HUVEC'
custom_gtf = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/HUVEC/HUVEC_merged_tama_gene_id_20_08_26.gtf'
isoquant_transcript_counts = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/HUVEC/HUVEC_quant/HUVEC/HUVEC.transcript_counts.tsv'


################################################################################
# HELPER FUNCTIONS
################################################################################

def count_junctions(bam_path, min_anchor=8, min_intron=40, min_mapq=1):
    """
    Get coordinates of junctions from bam files and count how often 
    they appear: get splice coordinates by parsing CIGAR string
    0-based half-open (chrom, start, end) -> read count.
    """
    counts = Counter()
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            # pos gives position on reference chromosome
            pos = read.reference_start
            # block used to anchor alignments: are they trustworthy (long enough at
            # both sites of the splice)
            block = 0
            blocks, introns = [], []
            for op, ln in read.cigartuples:
                if op in (0, 2, 7, 8):      # M, D, =, X
                    block += ln
                    pos += ln
                elif op == 1:               # I — insertion, consumes query only
                    continue
                elif op == 3:               # N — designates the intron bases
                    # alignment block to the left
                    blocks.append(block)
                    # write down genmoic coordinates of the junction: intron start, intron end
                    introns.append((pos, pos + ln))
                    pos += ln
                    # block to the right starst from 0
                    block = 0
                # S/H ignored: soft-clipped, hard-clipped, only present in read, not ref
            blocks.append(block)

            for i, (s, e) in enumerate(introns):
                if e - s < min_intron:
                    continue
                if blocks[i] >= min_anchor and blocks[i + 1] >= min_anchor:
                    counts[(read.reference_name, s, e)] += 1
    return counts


def junction_support(counts, all_junctions, tol=2):
    """
    For each GTF junction: (exact support, support within +/-tol).
    Returns dict: support[(chrom, s, e)]: (nr_jct_exact_support, nr_jct_support_sum_within_tol)
    """
    offsets = [(ds, de) for ds in range(-tol, tol + 1)
               for de in range(-tol, tol + 1)]
    support = {}
    for chrom, s, e in all_junctions:
        # sum all counts of BAM junctions within the tolerance window
        # get or return 0 if not present
        within = sum(counts.get((chrom, s + ds, e + de), 0)
                     for ds, de in offsets)
        # get the support for the exact junction too
        support[(chrom, s, e)] = (counts.get((chrom, s, e), 0), within)
    return support


def main(
    classification_txt,
    custom_gtf,
    bam_file,
    cell_type,
    outdir,
    isoquant_transcript_counts
):
    ############################################################################
    # READ IN DATA
    ############################################################################
    classification_df = pd.read_csv(
        classification_txt, sep='\t', header=0, index_col=None)

    isoquant_df = pd.read_csv(
        isoquant_transcript_counts, sep='\t', header=0, index_col=None)

    isoquant_dict = dict(zip(isoquant_df['feature_id'], isoquant_df['count']))

    classification_df['LR_isoquant'] = classification_df['isoform'].map(
        isoquant_dict)

    ############################################################################
    # GET LR SUPPORT
    ############################################################################

    # Get junction support of the long reads
    counts = count_junctions(bam_file)
    TID = re.compile(r'transcript_id "([^"]+)"')
    exons = defaultdict(list)
    chroms = {}

    # get junction coordinates in GTF file
    with open(custom_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            entry = line.rstrip("\n").split("\t")
            if entry[2] != "exon":
                continue
            transcript_match = TID.search(entry[8])
            if transcript_match is None:
                continue
            tx = transcript_match.group(1)
            exons[tx].append((int(entry[3]), int(entry[4]))
                             )   # 1-based inclusive
            chroms[tx] = entry[0]

    for tx, ex in exons.items():
        # sorting based on first coordinates (start), if these tie, then the second
        ex.sort()
        for i in range(len(ex) - 1):
            assert ex[i][1] < ex[i +
                                 1][0], f"{tx}: overlapping or duplicate exons {ex}"

    # GTF coordinates: (chr, exon end = intron start, exon start - 1 = intron end)
    # 1-based is GTF while BAM is 0-based half-open
    tx_junctions = {}          # transcript -> list of junction keys
    all_junctions = set()      # flat set of every annotated junction

    for tx, ex in exons.items():
        chrom = chroms[tx]
        # make a list of all junctions for this transcript, need to adjust to 0-based
        js = [(chrom, ex[i][1], ex[i + 1][0] - 1) for i in range(len(ex) - 1)]
        tx_junctions[tx] = js
        # update adds all junctions in list to set of GTF based junctions
        all_junctions.update(js)

    counts = count_junctions(bam_file)
    support = junction_support(counts, all_junctions, tol=2)

    min_cov = {}
    for tx, js in tx_junctions.items():
        if not js:
            min_cov[tx] = (None, None)              # mono-exonic
        else:
            min_cov[tx] = (min(support[j][0] for j in js),
                           min(support[j][1] for j in js))

    # get transcripts with minimum LR support
    transcript_with_lr_support = {t: counts[1] for t, counts in min_cov.items(
    ) if counts[1] is not None and counts[1] > 0}

    ############################################################################
    # FILTER CLASSIFICATION DF
    ############################################################################

    classification_df['min_LR_cov'] = classification_df['isoform'].map(
        transcript_with_lr_support)

    # SR supported transcripts:
    junction_supported_trans = classification_df[(classification_df["min_cov"] > 0) | (
        classification_df['min_LR_cov'] > 0)]

    mono_exon_supported_df = classification_df[((classification_df["subcategory"] == "mono-exon") |
                                                (classification_df["subcategory"] == "mono-exon_by_intron_retention")) & (classification_df['LR_isoquant'] > 2)]

    transcripts_filtered_df = pd.concat(
        [junction_supported_trans, mono_exon_supported_df], ignore_index=True)

    ############################################################################
    # WRITE OUTPUT FILES
    ############################################################################

    # write full classification_df in outdirectory
    # write transcripts_filtered_df
    # filter GTF and write


if __name__ == "__main__":
    args = parse_arguments()

    classification_txt = args.classification_txt
    custom_gtf = args.custom_gtf
    bam_file = args.bam_file
    cell_type = args.cell_type
    outdir = args.outdir
    isoquant_transcript_counts = args.isoquant_transcript_counts

    main(classification_txt, custom_gtf, bam_file, cell_type,
         outdir, isoquant_transcript_counts)
