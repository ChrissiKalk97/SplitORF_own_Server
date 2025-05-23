#!/bin/bash

#----- This script performs the preprocessing steps for the short RNA-seq samples ----- #
# ----- first cutadapt is used to trim adapters and filter for quality            ----- #
# ----- FASTQC is run to obtain quality control metrics and plots                 ----- #

eval "$(conda shell.bash hook)"
conda activate Riboseq

raw_data_dir="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData"
raw_data_fastqc_dir="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData/fastqc"

multiQC_outname="multiqc"
# 40 FASTQ files are indicated in the MD5 checksums

outidr_fastp="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/fastp"
STAR_out_dir="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/STAR"
STAR_index="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/STAR/index/Ens110"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
genome_annotation="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"


# bash fastqc_multiqc_RNA_seq.sh ${raw_data_dir} ${raw_data_fastqc_dir} ${multiQC_outname} raw

# bash trim_adapters_RNA_seq.sh ${raw_data_dir} ${outidr_fastp} _fastp

# bash fastqc_multiqc_RNA_seq.sh ${outidr_fastp} ${outidr_fastp}/fastqc ${multiQC_outname} fastp

bash STAR_map_short_RNA_seq.sh \
 $outidr_fastp \
 $STAR_out_dir \
 $STAR_index \
 $genome_fasta \
 $genome_annotation





