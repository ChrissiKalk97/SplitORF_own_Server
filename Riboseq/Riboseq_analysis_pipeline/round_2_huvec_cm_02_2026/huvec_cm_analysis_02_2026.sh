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
ensembl_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
bowtie_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"





huvec_dir="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_2_Feb_2026/uf_muellermcnicoll_2026_02_despic_RiboSeq_HUVEC_CM/HUVEC"
cm_dir="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_2_Feb_2026/uf_muellermcnicoll_2026_02_despic_RiboSeq_HUVEC_CM/CM"

out_dir="/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2"
out_dir_cm="${out_dir}"/CM
out_dir_huvec="${out_dir}"/HUVEC

report_dir="/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/round_2_huvec_cm_02_2026/outreports_of_runs"
report_dir_cm="${report_dir}"/CM
report_dir_huvec="${report_dir}"/HUVEC

echo "$data_dir"
echo "$out_dir"


######### Preprocessing: cutadapt, UMI trim, fastp #################################
# -d : dedup
# -c : cutadapt
# -p : paired end
# bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/Riboseq_pipeline.sh \
#  -d -p -u -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -i "$huvec_dir" -o "$out_dir_huvec" \
#  -r "$report_dir_huvec"

# bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/Riboseq_pipeline.sh \
#  -d -p -u -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -i "$cm_dir" -o "$out_dir_cm" \
#  -r "$report_dir_cm"

######### Align to transcriptome with bowtie1 #######################################
# bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/Riboseq_pipeline.sh \
#  -d -t "${bowtie_ref_fasta}" -i "$huvec_dir" -o "$out_dir_huvec" \
#  -r "$report_dir_huvec"

# bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/Riboseq_pipeline.sh \
#  -d -t "${bowtie_ref_fasta}" -i "$cm_dir" -o "$out_dir_cm" \
#  -r "$report_dir_cm"


######### Had wrong naming for deduplication bams, correct naming ####################
bash rename_dedup_bam_files.sh \
"${out_dir}/CM/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup"

bash rename_dedup_bam_files.sh \
"${out_dir}/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup"

conda activate r-env

Rscript /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/Ribowaltz/RiboWaltz_Michi_Vlado_2_Ingolia_HUVEC_CM_Feb2026.R \
"${out_dir}/CM/alignment_concat_transcriptome_Ignolia/Ribowaltz" \
"${out_dir}/CM/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup"

Rscript /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/Ribowaltz/RiboWaltz_Michi_Vlado_2_Ingolia_HUVEC_CM_Feb2026.R \
"${out_dir}/HUVEC/alignment_concat_transcriptome_Ignolia/Ribowaltz" \
"${out_dir}/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup"