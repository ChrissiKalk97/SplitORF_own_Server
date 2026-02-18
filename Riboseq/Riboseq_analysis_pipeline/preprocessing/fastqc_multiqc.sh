#!/bin/bash

INDIR=$1
OUTDIR_FASTQC=$2
multiQC_outname=$3


for FQ in "${INDIR}"/*.fastq.gz
do
    (fastqc \
    -o ${OUTDIR_FASTQC}/ \
    -t 32\
    ${FQ}) &
done

wait

multiqc --force --filename ${OUTDIR_FASTQC}/${multiQC_outname} ${OUTDIR_FASTQC}

