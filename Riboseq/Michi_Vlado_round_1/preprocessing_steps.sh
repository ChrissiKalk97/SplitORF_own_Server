#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

################################################################################
# RUN FASTQ                                                                    #
################################################################################
INDIR="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_1"
OUTDIR_FASTQC1="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/fastqc_unprocessed"

# ################################################################################
# # RUN FASTQC                                                                   #
# ################################################################################
# for FQ in "${INDIR}"/*.fastq.gz
# do
#     fastqc \
#     -o ${OUTDIR_FASTQC1}/ \
#     -t 32\
#     ${FQ}
# done


################################################################################
# RUN FASTP                                                                    #
################################################################################
OUTDIR_FASTP="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/fastp"

# for FQ in "${INDIR}"/*R1.fastq.gz
# do
# SAMPLE=$(basename "$FQ" .R1.fastq.gz)
# fastp \
#     --in1 ${FQ} \
#     --in2 ${INDIR}/${SAMPLE}.R2.fastq.gz \
#     --out1 ${OUTDIR_FASTP}/${SAMPLE}.R1.fastp.fastq.gz \
#     --out2 ${OUTDIR_FASTP}/${SAMPLE}.R2.fastp.fastq.gz \
#     --json ${OUTDIR_FASTP}/${SAMPLE}.fastp.json \
#     -h ${OUTDIR_FASTP}/${SAMPLE}.fastp.html \
#     --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#     --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     --overrepresentation_analysis \
#     --thread 32 \
#     --length_required 8 \
#     --correction \
#     --overlap_len_require 20

# UMI extraction is performed as the very first step, so the UMI things need to be run later
# -U \
# --umi_len=9 \
# --umi_loc=per_read
# done
# Parameters:
# default quality filter: Q15
# length required: not sure how long they will be after adapter
# trimming, adapters are also in the middle sometimes
# detect_adapter_for_pe: default would be overlap, not sure if there is overlap
# correction: PE reads correct overlap (if there is) mismathc by the base with higher quality
# U: enable UMi processing, first 9 bp, append UMI1_UMI2 to the read names and trim both reads of the sequences
# enable analysis of overrepresented sequences!


# ################################################################################
# # RUN FASTQC AFTER FASTP                                                       #
# ################################################################################

fastp_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/fastp"

for FQ in "${fastp_dir}"/*.fastq.gz
do
    fastqc \
    -o ${OUTDIR_FASTP}/fastqc/ \
    -t 32\
    ${FQ}
done


