#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

# ensembl_filtered_gtf=$1
# genome_fasta=$2
# consensus_reads_fofn=$3


reference_gtf="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
long_read_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"
script_dir="/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3"
out_path="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3"


cell_type=HUVEC
consensus_reads_fofn="pacbio_consensus_${cell_type}.fofn"
bam_dir="/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment/${cell_type}/pbmm2_align"
short_read_dir="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/${cell_type}_fastp"


bash stringtie3_correct_ref_filtered_LRs_05_01_26.sh \
    -b "${bam_dir}" \
    -c "$cell_type" \
    -e "${ensembl_full_gtf}" \
    -f "${genome_fasta}" \
    -l "${long_read_dir}" \
    -o "${out_path}" \
    -p "${short_read_dir}" \
    -n "${consensus_reads_fofn}" \
    -r "${reference_gtf}" \
    -s "${script_dir}" 



cell_type=CM
consensus_reads_fofn="pacbio_consensus_${cell_type}.fofn"
bam_dir="/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment/${cell_type}/pbmm2_align"
short_read_dir="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/${cell_type}_fastp"

bash stringtie3_correct_ref_filtered_LRs_05_01_26.sh \
    -b "${bam_dir}" \
    -c "$cell_type" \
    -e "${ensembl_full_gtf}" \
    -f "${genome_fasta}" \
    -l "${long_read_dir}" \
    -o "${out_path}" \
    -p "${short_read_dir}" \
    -n "${consensus_reads_fofn}" \
    -r "${reference_gtf}" \
    -s "${script_dir}" 