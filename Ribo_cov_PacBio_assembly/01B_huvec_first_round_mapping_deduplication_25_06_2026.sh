#!/bin/bash


# =============================================================================
# Script Name: run_Riboseq_analysis.sh
# Description: This script performs a custom analysis pipeline for the Riboseq 
#               analysis of data created with Vlado's protocol (April 2025):
#              - Step 1: Data preprocessing
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
module_dir="/home/ckalk/scripts/Ribo_seq_analysis/Riboseq_analysis_pipeline"

genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_23_June_2026/HUVEC/HUVEC_merged_tama_gene_id.gtf"
star_index_huvec="/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2_PacBio_assembly/HUVEC/alignment_genome/STAR/index"

indir="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_1"
OUTDIR_FASTQC1="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/fastqc_unprocessed"
OUTDIR_CUTADAPT="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt"
Umi_adpt_trimmed_path="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/UMI_trimmed_custom"
FASTP_OUT="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim"
fastpFASTQC="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim/fastqc"

out_dir="/projects/splitorfs/work/Riboseq/Output/HUVEC_round_1_PacBio_assembly"
out_dir_huvec="${out_dir}"/HUVEC

mkdir -p "$out_dir_huvec"

report_dir="/home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/outreports_of_runs"
report_dir_huvec="${report_dir}"/HUVEC

mkdir -p "$report_dir_huvec"



######### Align to genome with STAR and deduplicate #############################
# start with the preprocessed reads:
# as done in here for the previous analysis:
# /home/ckalk/scripts/Ribo_seq_analysis/Riboseq_analysis_pipeline/run_Riboseq_analysis.sh
 

bash "${module_dir}"/Riboseq_pipeline_2.sh \
 -a "${star_index_huvec}" -i "$FASTP_OUT" -o "$out_dir_huvec" -g $huvec_gtf -m "${module_dir}" -d -s -p -q \
  > "$report_dir_huvec"/Riboseq_pipeline_mapping_first_round_HUVEC_Mando_Iso_Stringtie_rescue_Assembly_24_06_2026.out 2>&1

