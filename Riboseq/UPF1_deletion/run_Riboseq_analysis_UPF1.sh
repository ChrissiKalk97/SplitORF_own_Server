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

data_dir="/projects/splitorfs/work/UPF1_deletion/data"
out_dir="/projects/splitorfs/work/UPF1_deletion/Output"
module_dir="/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_analysis_pipeline"


Bowtie2_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
Bowtie2_base_name="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/index"
Bowtie2_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"


UMI_Indir_transcriptomic="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"
UMI_dedup_outdir_transcriptomic="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/deduplicated"

echo "$data_dir"
echo "$out_dir"

# get FASTQ files
# bash "${module_dir}"/Riboseq_pipeline.sh -i "$data_dir" -o "$out_dir" -m "${module_dir}" -b

data_dir="${data_dir}/fastq"
bash "${module_dir}"/Riboseq_pipeline.sh \
 -i "$data_dir" -o "$out_dir" -g $ensembl_gtf -m "${module_dir}" -q -d