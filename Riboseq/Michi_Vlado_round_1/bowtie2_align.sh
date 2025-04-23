#!/bin/bash

Bowtie2_base_name=$1
Bowtie_ref_fasta=$2
FASTQ_Inpath=$3
Bowtie_out_dir=$4
aligned_name=$5

# bowtie2-build ${Bowtie_ref_fasta} ${Bowtie2_base_name} --threads 32


for FQ in "${FASTQ_Inpath}"/*R1*.fastq.gz
do
    sample=$(basename "$FQ")       # remove path
    sample=${sample%%R1*}          # remove R1 and everything after
    echo ${sample}
    FQ2=${FQ/R1/R2} # substitute R1 with R2
    bowtie2 \
    -p 32 \
    -L 20 \
    --very-sensitive-local \
    -N 1 \
    --no-discordant \
    -S "${Bowtie_out_dir}"/"${sample}"bowtie2_${aligned_name}.sam \
    --un-conc-gz "${Bowtie_out_dir}"/"${sample}"bowtie2_tRNA_unaligned\
    -x ${Bowtie2_base_name} \
    -q \
    -1 ${FQ} \
    -2 ${FQ2}
done