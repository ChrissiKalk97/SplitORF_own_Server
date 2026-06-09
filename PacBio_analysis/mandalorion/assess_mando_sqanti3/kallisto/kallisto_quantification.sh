#!/bin/bash


eval "$(conda shell.bash hook)"
conda activate Riboseq

kallisto_index=$1
fastp_indir=$2
kallisto_outdir=$3

echo "kallisto params"
echo $kallisto_index
echo $fastp_indir
echo $kallisto_outdir

for reads in ${fastp_indir}/*R1.fastp.fastq.gz; do

     echo $reads

     sample=$(basename "$reads" _merged_fastp.R1.fastp.fastq.gz)

     echo $sample

     kallisto quant -i ${kallisto_index}\
     -o $kallisto_outdir/${sample}\
     -b 100 \
     -t 16 \
     $reads \
     ${fastp_indir}/${sample}_merged_fastp.R2.fastp.fastq.gz
done
