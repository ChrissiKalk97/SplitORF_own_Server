#!/bin/bash

#----- This script performs the preprocessing steps for the short RNA-seq samples ----- #
# ----- first cutadapt is used to trim adapters and filter for quality            ----- #
# ----- FASTQC is run to obtain quality control metrics and plots                 ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio

raw_data_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData"
merged_data_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/merged"
raw_data_fastqc_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData/fastqc"

multiQC_outname="multiqc"
# 40 FASTQ files are indicated in the MD5 checksums

outdir_fastp="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_fastp"
Genome_Fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
EnsemblFilteredRef="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
SalmonOutDir="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Salmon/Ens_110"
Decoys="/projects/splitorfs/work/reference_files/decoys.txt"

# bash analyze_short_RNA_seq_cell_type.sh \
#     "${raw_data_dir_huvec}" \
#     "${merged_data_dir_huvec}" \
#     "${raw_data_fastqc_dir_huvec}" \
#     "${outdir_fastp}"

# bash Salmon_quantification.sh \
#   ${Genome_Fasta} \
#   ${EnsemblFilteredRef} \
#   ${SalmonOutDir} \
#   ${outdir_fastp} \
#   ${SalmonOutDir} \
#   ${Decoys}

conda activate r-env
Rscript DeSeq2_Salmon_counts_HUVEC_Ens_110_for_Ranjan.R

