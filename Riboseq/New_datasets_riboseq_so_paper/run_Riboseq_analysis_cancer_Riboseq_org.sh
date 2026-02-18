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

module_dir="/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_analysis_pipeline"



meta_stat_breast_cancer_data_dir="/projects/splitorfs/work/Riboseq/data/cancer_data_riboseq_org/metastat_breast_cancer_PRJNA898352/fastq"
glioblastoma_data_dir="/projects/splitorfs/work/Riboseq/data/cancer_data_riboseq_org/glioblastoma_PRJNA591767/fastq"
breast_cancer_data_dir="/projects/splitorfs/work/Riboseq/data/cancer_data_riboseq_org/breast_cancer_PRJNA523167/fastq"

heart_dir="/projects/splitorfs/work/Riboseq/data/heart_iPSC/fastq"
leukemia_dir="/projects/splitorfs/work/Riboseq/data/leukemia/fastq"
endothel_dir="/projects/splitorfs/work/Riboseq/data/endothel_Siragusa/fastq_from_bam"

out_dir="/projects/splitorfs/work/Riboseq/data"

echo "$data_dir"
echo "$out_dir"


# -c cutadapt for general preprocessing, misleading does use fastp...
# do not run the metastatic breast cancer because the data is weird, just exclude from analysis
# bash "${module_dir}"/Riboseq_pipeline.sh \
# -i "$meta_stat_breast_cancer_data_dir" -o "$out_dir" \
# -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -f 5 -m "${module_dir}"

# bash "${module_dir}"/Riboseq_pipeline.sh -i "$glioblastoma_data_dir" \
# -o "$out_dir" -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m "${module_dir}"

bash "${module_dir}"/Riboseq_pipeline.sh -i "$breast_cancer_data_dir" \
 -o "$out_dir" -c CTGTAGGCACCATCAAT -m "${module_dir}"

# bash "${module_dir}"/Riboseq_pipeline.sh -i "$heart_dir" \
# -o "$out_dir" -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m "${module_dir}"

# bash "${module_dir}"/Riboseq_pipeline.sh -i "$leukemia_dir" \
# -o "$out_dir" -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m "${module_dir}"

# bash "${module_dir}"/Riboseq_pipeline.sh -i "$endothel_dir" \
# -o "$out_dir" -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m "${module_dir}"