"""Filter a GTF by LR and SR support (after SQANTI QC).

Counts splice junctions in a long-read BAM, matches them to the junctions
implied by the custom GTF (allowing JUNCTION_TOL bp of wobble), and keeps
transcripts with either short-read or long-read junction support.
Mono-exonic transcripts have no junctions and are kept on IsoQuant
expression instead.
This is performed for the TAMA merged assemblies of Isoquant, Strigntie3 and 
Mandalorion after the SQANTI pipeline of QC, Filter and Rescue.
Usage:
    python filter_by_lr_support.py \\
        --classification_txt classification.txt (from SQANTI QC) \\
        --custom_gtf custom.gtf \\
        --bam_file merged.bam \\
        --isoquant_transcript_counts transcript_counts.tsv \\
        --celltyppe HUVEC or CM \\
        --outdir results/

Outputs:
    <outdir>/<cell_type>_filtered_classification.txt
    <outdir>/<cell_type>_filtered.gtf
"""


import os
import pandas as pd
import argparse
import pysam
from collections import Counter
from collections import defaultdict
import re
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pkg_resources")

# for the isoquant conda environment


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process peptide and SO ID mapping files.")

    parser.add_argument(
        "--classification_txt",
        required=True,
        help="Path to SQANTI QC classification file."
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
# COSTANTS
################################################################################
# GTF attribute parsing
TID_RE = re.compile(r'(?:^|;)\s*transcript_id "([^"]+)"')
GID_RE = re.compile(r'(?:^|;)\s*gene_id "([^"]+)"')

# Filter thresholds
JUNCTION_TOL = 2            # bp of wobble allowed matching BAM to GTF junctions
MIN_JUNCTION_READS = 1      # LR reads at the least-covered junction
MIN_MONO_EXON_READS = 3     # IsoQuant counts for mono-exonic transcripts


################################################################################
# PATH DEFINITIONS FOR TESTING
################################################################################
# classification_txt = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/SQANTI3_QC/HUVEC/isoforms_classification.txt'
# cell_type = 'HUVEC'
# bam_file = '/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment/HUVEC/minimap2_align/merged/HUVEC_merged_sorted.bam'
# outdir = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/HUVEC'
# custom_gtf = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/HUVEC/HUVEC_merged_tama_gene_id_20_08_26.gtf'
# isoquant_transcript_counts = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/HUVEC/HUVEC_quant/HUVEC/HUVEC.transcript_counts.tsv'


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


def parse_gtf_exons(gtf_path):
    """Extract per-transcript exon coordinates from a GTF.

    Exons are sorted by start position. Coordinates are 1-based inclusive,
    as in the GTF itself — conversion to 0-based happens in build_junctions.
    Get exon start end coordinates per transcript as dict and
    dict of which transcripts belong to which chromosomes

    Raises:
        ValueError: if a transcript has overlapping or duplicate exons.
    Returns:
        transcript_exon_dict (exons), transcript_chromosome_dict (chroms)
    """
    exons = defaultdict(list)
    chroms = {}

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, _, attr = line.rstrip(
                "\n").split("\t")
            if feature != "exon":
                continue
            transcript_match = TID_RE.search(attr)
            if transcript_match is None:
                continue
            tx = transcript_match.group(1)
            exons[tx].append((int(start), int(end)))   # 1-based inclusive
            chroms[tx] = chrom

    for tx, ex in exons.items():
        # sorting based on first coordinates (start), if these tie, then the second
        ex.sort()
        for i in range(len(ex) - 1):
            if ex[i][1] >= ex[i + 1][0]:
                raise ValueError(f"{tx}: overlapping or duplicate exons {ex}")
    return exons, chroms


def build_junctions(exons, chroms):
    """Convert per-transcript exon lists into junction coordinates.

    GTF exons are 1-based inclusive; junctions are returned 0-based
    half-open to match pysam, so a junction is
    (chrom, exon_i.end, exon_i+1.start - 1).

    Returns:
        tx_junctions: transcript -> list of junction keys
        all_junctions: flat set of every annotated junction
    """
    tx_junctions = {}  # transcript -> list of junction keys
    all_junctions = set()  # flat set of every annotated junction
    for tx, ex in exons.items():
        chrom = chroms[tx]
        # make a list of all junctions for this transcript, need to adjust to 0-based
        js = [(chrom, ex[i][1], ex[i + 1][0] - 1) for i in range(len(ex) - 1)]
        tx_junctions[tx] = js
        # update adds all junctions in list to set of GTF based junctions
        all_junctions.update(js)
    return tx_junctions, all_junctions


def min_junction_coverage(tx_junctions, support):
    """Least-covered junction per transcript, as (exact, within-tolerance).

    Mono-exonic transcripts get (None, None) — they have no junctions, so
    the test does not apply to them.
    """
    min_cov = {}
    for tx, js in tx_junctions.items():
        if not js:
            min_cov[tx] = (None, None)
        else:
            min_cov[tx] = (min(support[j][0] for j in js),
                           min(support[j][1] for j in js))
    return min_cov


def filter_gtf(gtf_in, gtf_out, keep_tx):
    """Keep lines whose transcript is in keep_tx; keep genes that retain >=1."""
    keep_genes = set()
    all_tx = set()
    all_genes = set()
    with open(gtf_in) as fh:                # pass 1, get genes to keep
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] == "gene":
                g = GID_RE.search(f[8])
                if g:
                    all_genes.add(g.group(1))
                continue
            m = TID_RE.search(f[8])
            if m is None:
                continue
            all_tx.add(m.group(1))
            if m and m.group(1) in keep_tx:
                g = GID_RE.search(f[8])
                if g:
                    keep_genes.add(g.group(1))

    n_out = 0
    with open(gtf_in) as fh, open(gtf_out, "w") as out:   # pass 2: filter and write
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.split("\t")
            if f[2] == "gene":
                g = GID_RE.search(f[8])
                if g and g.group(1) in keep_genes:
                    out.write(line)
                    n_out += 1
            else:
                m = TID_RE.search(f[8])
                if m and m.group(1) in keep_tx:
                    out.write(line)
                    n_out += 1

    print(f"genes: {len(all_genes)} -> {len(keep_genes)}")
    print(f"transcripts: {len(all_tx)} -> {len(keep_tx & all_tx)}")
    missing = keep_tx - all_tx
    # need to check in case classification TXT and GTF differ in their transcripts
    if missing:
        print(f"warning: {len(missing)} kept transcripts absent from GTF, "
              f"e.g. {sorted(missing)[:3]}")
    return n_out


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
        classification_txt, sep='\t', header=0, index_col=None, dtype={'chrom': str})

    isoquant_df = pd.read_csv(
        isoquant_transcript_counts, sep='\t', header=0, index_col=None)

    isoquant_dict = dict(zip(isoquant_df['feature_id'], isoquant_df['count']))

    classification_df['LR_isoquant'] = classification_df['isoform'].map(
        isoquant_dict)

    ############################################################################
    # GET LR SUPPORT
    ############################################################################
    # get exon dict for transcripts (coordinates), get chrom dict for transcripts
    exons, chroms = parse_gtf_exons(custom_gtf)

    tx_junctions, all_junctions = build_junctions(exons, chroms)

    counts = count_junctions(bam_file)
    support = junction_support(counts, all_junctions, tol=2)

    min_cov = min_junction_coverage(tx_junctions, support)

    # get transcripts with minimum LR support
    lr_cov = {t: cov[1] for t, cov in min_cov.items() if cov[1] is not None}
    # 0 means "no support", NaN means "mono-exonic or not in GTF".
    classification_df['min_LR_cov'] = classification_df['isoform'].map(lr_cov)

    # add exact junction coverage to classification TXT
    lr_cov_exact = {t: cov[0]
                    for t, cov in min_cov.items() if cov[0] is not None}
    classification_df['min_LR_cov_exact'] = classification_df['isoform'].map(
        lr_cov_exact)

    ############################################################################
    # FILTER CLASSIFICATION DF
    ############################################################################
    mono = classification_df["subcategory"].isin(
        ["mono-exon", "mono-exon_by_intron_retention"])

    keep = (
        (classification_df["min_cov"] >= MIN_JUNCTION_READS)
        | (classification_df["min_LR_cov"] >= MIN_JUNCTION_READS)
        | (mono & (classification_df["LR_isoquant"] >= MIN_MONO_EXON_READS))
    )
    transcripts_filtered_df = classification_df[keep].copy()

    assert not transcripts_filtered_df['isoform'].duplicated().any()

    # additional checks
    print(
        f"{classification_df['LR_isoquant'].isna().sum()} transcripts absent from IsoQuant")
    # if SQANTI run w/o short reads, then this script does not make sense
    if classification_df['min_cov'].isna().all():
        raise ValueError(
            "min_cov is entirely NA — SQANTI3 was run without short-read coverage, "
            "so the SR branch of the filter cannot fire"
        )

    ############################################################################
    # WRITE OUTPUT FILES
    ############################################################################
    # write full classification_df in outdirectory
    classification_df.to_csv(os.path.join(
        outdir, f'{cell_type}_LR_SR_support_classification.txt'), sep='\t', index=False)
    # write transcripts_filtered_df
    transcripts_filtered_df.to_csv(os.path.join(
        outdir, f'{cell_type}_LR_SR_support_classification_filtered.txt'), sep='\t', index=False)
    # filter GTF and write
    keep_tx = set(transcripts_filtered_df['isoform'])
    filter_gtf(custom_gtf,
               os.path.join(
                   outdir, f'{cell_type}_LR_SR_support_filtered.gtf'),
               keep_tx)


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
