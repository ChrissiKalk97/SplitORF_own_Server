#!/bin/bash

indir=$1
outdir_fastp=$2
append_outfile=$3


for FQ in ${indir}/*_1.fq.gz; 
do
sample=$(basename "$FQ")
sample=${sample%%_1*}     
FQ2=${FQ/_1.fq.gz/_2.fq.gz}
fastp \
    --in1 ${FQ} \
    --in2 ${FQ2} \
    --out1 ${outdir_fastp}/${sample}${append_outfile}.R1.fastp.fastq.gz \
    --out2 ${outdir_fastp}/${sample}${append_outfile}.R2.fastp.fastq.gz \
    --json ${outdir_fastp}/${sample}fastp.json \
    -h ${outdir_fastp}/${sample}fastp.html \
    --adapter_sequence  AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --trim_poly_g \
    --trim_poly_x \
    --overrepresentation_analysis \
    --thread 16 \
    --length_required 25 \
    --correction \
    --overlap_len_require 30
done

multiqc --force --filename ${outdir_fastp}/fastp_multiqc_report ${outdir_fastp}