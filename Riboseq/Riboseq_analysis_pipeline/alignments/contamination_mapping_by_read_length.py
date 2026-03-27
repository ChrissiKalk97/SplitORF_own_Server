import pysam
import sys
import os
import os.path
from Bio import SeqIO
import pandas as pd


def assign_trans_IDs_to_source(ref_fastas):
    """
    Assign the IDs in reference fasta files to the source

    Return a dict with the sources as keys and the transcript IDs as values

    """
    ref_fastas = ref_fastas.split(",")
    reference_dict = {}
    for fasta in ref_fastas:
        ref_seqs = SeqIO.parse(fasta, "fasta")
        source = os.path.basename(fasta)
        source = source.rsplit('.', 1)[0]
        for seq in ref_seqs:
            reference_dict[seq.id] = source
    return reference_dict


ref_fastas = sys.argv[1]
bam_dir = sys.argv[2]

# ref_fastas = '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta'
# bam_dir = '/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/filtered'

bam_files = []

for filename in os.listdir(bam_dir):
    if filename.endswith(".bam"):
        bam_files.append(os.path.join(bam_dir, filename))

out_dir = os.path.dirname(bam_files[0])
os.mkdir(os.path.join(out_dir, 'read_mapping_by_length'))

reference_dict = assign_trans_IDs_to_source(ref_fastas)
summary_mapping_stats_df = None

results = []
for bam_file in bam_files:
    bam = pysam.AlignmentFile(bam_file)
    bam_name = os.path.basename(bam_file).split('.')[0]
    total_reads = 0
    count_dict = {}

    for read in bam:
        # skip unwanted reads
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        total_reads += 1
        length = read.query_length
        transcript = read.reference_name

        if transcript is None:
            continue

        # default to mRNA if missing
        source = reference_dict[transcript]
        if length in count_dict.keys():
            if source in count_dict[length].keys():
                count_dict[length][source] += 1
            else:
                count_dict[length][source] = 1
        else:
            count_dict[length] = {}
            count_dict[length][source] = 1

    mapping_df = pd.DataFrame(count_dict).fillna(0).astype(int)
    mapping_df = mapping_df.transpose()
    mapping_df = mapping_df.sort_index()
    mapping_df['total'] = mapping_df.sum(axis=1)
    mapping_df_percent = mapping_df.copy()
    for col in mapping_df.columns:
        mapping_df_percent[f'{col} %'] = mapping_df[col]/mapping_df['total']
    mapping_df_percent.index.name = "length"
    mapping_df_percent = mapping_df_percent.reset_index()
    mapping_df_percent.to_csv(os.path.join(
        out_dir, 'read_mapping_by_length', f'{bam_name}_mapping_by_read_length.csv'))
