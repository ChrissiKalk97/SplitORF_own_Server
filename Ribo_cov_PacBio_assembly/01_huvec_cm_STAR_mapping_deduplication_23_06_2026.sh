#!/bin/bash


# =============================================================================
# Script Name: run_Riboseq_analysis.sh
# Description: This script performs a custom analysis pipeline for the Riboseq 
#               analysis of data created with Vlado's protocol (April 2025):
#              - Step 1: Data preprocessing (skipped)
#              - Step 2: Transcriptomic alignment (Ingolia reference)
#              - Step 3: Transcriptomic deduplication
#              - Step 4: Genomic alignment, deduplication
#              - Step 5: Split-ORF analysis
# Usage:       bash my_pipeline.sh 
# Author:      Christina Kalk
# Date:        2025-05-27
# =============================================================================






eval "$(conda shell.bash hook)"
conda activate Riboseq

################################################################################
# PATH DEFINTIONS                                                              #
################################################################################

genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_23_June_2026/HUVEC/HUVEC_merged_tama_gene_id.gtf"
cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_23_June_2026/CM/CM_merged_tama_gene_id.gtf"
star_index_cm="/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2_PacBio_assembly/CM/alignment_genome/STAR/index"
star_index_huvec="/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2_PacBio_assembly/HUVEC/alignment_genome/STAR/index"


# update to the new location of the Ribo-seq preprocessing pipeline
module_dir="/home/ckalk/scripts/Ribo_seq_analysis/Riboseq_analysis_pipeline"


# huvec_dir="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_2_Feb_2026/uf_muellermcnicoll_2026_02_despic_RiboSeq_HUVEC_CM/HUVEC"
# cm_dir="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_2_Feb_2026/uf_muellermcnicoll_2026_02_despic_RiboSeq_HUVEC_CM/CM"

out_dir="/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2_PacBio_assembly"
out_dir_cm="${out_dir}"/CM
out_dir_huvec="${out_dir}"/HUVEC

mkdir -p "$out_dir_cm"
mkdir -p "$out_dir_huvec"

fastp_huvec_dir="/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/preprocess/cutadapt/fastp_filter_after_UMI_trim"
fastp_cm_dir="/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/CM/preprocess/cutadapt/fastp_filter_after_UMI_trim"

report_dir="/home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/outreports_of_runs"
report_dir_cm="${report_dir}"/CM
report_dir_huvec="${report_dir}"/HUVEC

mkdir -p "$report_dir_cm"
mkdir -p "$report_dir_huvec"

echo "$data_dir"
echo "$out_dir"


# These preprocessing steps were already performed for the initial analysis
# otherwise would have needed to perform them here
######### Preprocessing: cutadapt, UMI trim, fastp #################################
# -d : dedup
# -c : cutadapt
# -p : paired end
# bash "${module_dir}"/Riboseq_pipeline.sh \
#  -d -p -u -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -i "$huvec_dir" -o "$out_dir_huvec" \
#  -r "$report_dir_huvec" -m "${module_dir}"

# bash "${module_dir}"/Riboseq_pipeline.sh \
#  -d -p -u -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -i "$cm_dir" -o "$out_dir_cm" \
#  -r "$report_dir_cm" -m "${module_dir}"




######### Align to genome with STAR and deduplicate #############################
# need to put explicitly the fastp directory, do not need to do this for transcriptomic
# alignment
 
 #-q 
bash "${module_dir}"/Riboseq_pipeline.sh \
 -a "${star_index_huvec}" -i "$fastp_huvec_dir" -o "$out_dir_huvec" -g $huvec_gtf -m "${module_dir}" -d \
  > "$report_dir_huvec"/Riboseq_pipeline_HUVEC_Mando_Iso_Stringtie_rescue_Assembly_23_06_2026.out 2>&1


  bash "${module_dir}"/Riboseq_pipeline.sh \
 -a "${star_index_cm}" -i "$fastp_cm_dir" -o "$out_dir_cm" -g $cm_gtf -m "${module_dir}" -d \
 > "$report_dir_cm"/Riboseq_pipeline_deduplication_CM_Mando_Iso_Stringtie_rescue_Assembly_23_06_2026.out 2>&1





