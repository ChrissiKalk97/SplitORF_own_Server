#!/bin/bash

#----- This script performs the preprocessing steps for the short RNA-seq samples ----- #
# ----- first cutadapt is used to trim adapters and filter for quality            ----- #
# ----- FASTQC is run to obtain quality control metrics and plots                 ----- #

eval "$(conda shell.bash hook)"
conda activate Riboseq

raw_data_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData"
merged_data_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/merged"
raw_data_fastqc_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData/fastqc"

multiQC_outname="multiqc"
# 40 FASTQ files are indicated in the MD5 checksums

outidr_fastp="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_fastp"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
decoys="/projects/splitorfs/work/reference_files/decoys.txt"

out_dir_huvec="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_TAMA_08_10_26"
star_out_dir_huvec="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_TAMA_08_10_26/STAR"
salmon_out_dir_huvec_tama="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_TAMA_08_10_26/Salmon"
huvec_tama_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/HUVEC/HUVEC_LR_SR_support_filtered.gtf"

mkdir -p "/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Salmon"
mkdir -p "$out_dir_huvec"

bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/analyze_short_RNA_seq_cell_type.sh \
 "$raw_data_dir_huvec" \
 "$merged_data_dir_huvec" \
 "$raw_data_fastqc_dir_huvec" \
 "$outidr_fastp_huvec"

# bash STAR_map_short_RNA_seq.sh \
#  "$outidr_fastp" \
#  "$star_out_dir_huvec" \
#  "$star_out_dir_huvec/index" \
#  "$genome_fasta" \
#  "${huvec_tama_gtf}"


python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
    "$star_out_dir_huvec" \
    STAR_align_genome.csv


# bash Salmon_quantification.sh \
#   "${genome_fasta}" \
#   "${huvec_tama_gtf}" \
#   "${salmon_out_dir_huvec_tama}" \
#   "${outidr_fastp}" \
#   "${salmon_out_dir_huvec_tama}" \
#   "${decoys}"



#################################################################################
# ------------------ ANALYZE CM SHORT-READ DATA              ------------------ #
#################################################################################
raw_data_dir_cm="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/CM_short_reads/X208SC25032333-Z01-F002/01.RawData"
merged_data_dir_cm="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/CM_short_reads/merged"
raw_data_fastqc_dir_cm="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/CM_short_reads/X208SC25032333-Z01-F002/01.RawData/fastqc"

multiQC_outname="multiqc"
# 40 FASTQ files are indicated in the MD5 checksums

outidr_fastp_cm="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/CM_fastp"
out_dir_cm="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/CM_TAMA_08_10_26"
star_out_dir_cm="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/CM_TAMA_08_10_26/STAR"
salmon_out_dir_cm_tama="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/CM_TAMA_08_10_26/Salmon"


cm_tama_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/CM/CM_LR_SR_support_filtered.gtf"


mkdir -p "$out_dir_cm"

bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/analyze_short_RNA_seq_cell_type.sh \
 "$raw_data_dir_cm" \
 "$merged_data_dir_cm" \
 "$raw_data_fastqc_dir_cm" \
 "$outidr_fastp_cm"

# bash STAR_map_short_RNA_seq.sh \
#  "$outidr_fastp_cm" \
#  "$star_out_dir_cm" \
#  "$star_out_dir_cm/index" \
#  "$genome_fasta" \
#  "${cm_tama_gtf}"


python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
   "$star_out_dir_cm" \
    STAR_align_genome.csv


bash Salmon_quantification.sh \
  "${genome_fasta}" \
  "${cm_tama_gtf}" \
  "${salmon_out_dir_cm_tama}" \
  "${outidr_fastp_cm}" \
  "${salmon_out_dir_cm_tama}" \
  "${decoys}"


#################################################################################
# ------------------ PERFORM DESeq2 DEG ANALYSIS             ------------------ #
#################################################################################
conda activate r-env

Rscript /home/ckalk/scripts/SplitORFs/short_RNA_seq/DeSeq2_Salmon_Giibs_HUVEC_CM_12_10_26.R \
"${cm_tama_gtf}" \
"${salmon_out_dir_cm_tama}" \
"5NMD,DHYPO,DNOR" \
"1,2,3"


Rscript /home/ckalk/scripts/SplitORFs/short_RNA_seq/DeSeq2_Salmon_Giibs_HUVEC_CM_12_10_26.R \
"${huvec_tama_gtf}" \
"${salmon_out_dir_huvec_tama}" \
"50NMD,5NMD,DHYPO,DMSO,DNOR" \
"2,3,4"
