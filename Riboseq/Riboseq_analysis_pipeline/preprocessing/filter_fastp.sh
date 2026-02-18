#!/bin/bash

INDIR=$1
OUTDIR_FASTP=$2
append_outfile=$3


for FQ in "${INDIR}"/*R1.*fastq.gz
do
SAMPLE=$(basename "$FQ")
SAMPLE=${SAMPLE%%R1*}     
FQ2=${FQ/R1/R2}
fastp \
    --in1 ${FQ} \
    --in2 ${FQ2} \
    --out1 ${OUTDIR_FASTP}/${SAMPLE}${append_outfile}.R1.fastp.fastq.gz \
    --out2 ${OUTDIR_FASTP}/${SAMPLE}${append_outfile}.R2.fastp.fastq.gz \
    --json ${OUTDIR_FASTP}/${SAMPLE}fastp.json \
    -h ${OUTDIR_FASTP}/${SAMPLE}fastp.html \
    --disable_adapter_trimming \
    --overrepresentation_analysis \
    --thread 16 \
    --length_required 25 \
    --length_limit 50 \
    --correction \
    --overlap_len_require 20
done
# Parameters:
# default quality filter: Q15
# length required: not sure how long they will be after adapter
# trimming, adapters are also in the middle sometimes
# detect_adapter_for_pe: default would be overlap, not sure if there is overlap
# correction: PE reads correct overlap (if there is) mismathc by the base with higher quality
# U: enable UMi processing, first 9 bp, append UMI1_UMI2 to the read names and trim both reads of the sequences
# enable analysis of overrepresented sequences!