#!/bin/bash

#----- This script intersects the ultraconserved regions from literature with  ----- #
# ----- Split-ORf unique regions, stop codons            ----- #
# ----- FASTQC is run to obtain quality control metrics and plots                 ----- #

eval "$(conda shell.bash hook)"
conda activate Riboseq

raw_data_dir=$1
merged_data_dir=$2
raw_data_fastqc_dir=$3
outidr_fastp=$4

multiQC_outname="multiqc"

echo $raw_data_dir
echo $merged_data_dir
echo $raw_data_fastqc_dir
echo $outidr_fastp

if [ ! -e ${raw_data_fastqc_dir} ]; then
    bash fastqc_multiqc_RNA_seq.sh ${raw_data_dir} ${raw_data_fastqc_dir} ${multiQC_outname} raw
else
    echo ${raw_data_fastqc_dir} found skipping FASTQC step
fi

if [ ! -e ${merged_data_dir} ]; then
    bash merge_fastq_files.sh ${raw_data_dir} ${merged_data_dir}
else
    echo ${merged_data_dir} found skipping merging of short read files
fi

if [ ! -e ${merged_data_dir}/fastqc ]; then
    bash fastqc_multiqc_RNA_seq.sh ${merged_data_dir} ${merged_data_dir}/fastqc ${multiQC_outname} raw
else
    echo ${merged_data_dir}/fastqc found skipping FASTQC
fi

if [ ! -e ${outidr_fastp} ]; then
    mkdir ${outidr_fastp}
    bash trim_adapters_RNA_seq.sh ${merged_data_dir} ${outidr_fastp} _fastp
else
    echo ${outidr_fastp} found skipping fastp
fi

if [ ! -e ${outidr_fastp}/fastqc ]; then
    bash fastqc_multiqc_RNA_seq.sh ${outidr_fastp} ${outidr_fastp}/fastqc ${multiQC_outname} fastp
else
    echo ${outidr_fastp}/fastqc found skipping FASTQC
fi

