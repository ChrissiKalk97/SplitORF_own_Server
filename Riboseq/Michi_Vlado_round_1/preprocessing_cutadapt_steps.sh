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
# RUN CUTADAPT TO TRIM ADAPTERS                                                #
################################################################################
OUTDIR_CUTADAPT="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt"

# for FQ in "${INDIR}"/*R1.fastq.gz
# do
# SAMPLE=$(basename "$FQ" .R1.fastq.gz)

# cutadapt \
#     -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#     -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     -o ${OUTDIR_CUTADAPT}/${SAMPLE}.R1.cutadapt.fastq.gz \
#     -p ${OUTDIR_CUTADAPT}/${SAMPLE}.R2.cutadapt.fastq.gz \
#     --cores 16 \
#     --minimum-length 20 \
#     --max-n 0.1 \
#     --max-expected-errors 1 \
#     ${FQ} \
#     ${INDIR}/${SAMPLE}.R2.fastq.gz


   

# done
# Parameters:
# trim the adapter no matter where they are located and if they are parital 
# only trim at 3' end
# Require min length of 20 (inlcuding the UMI, could even be stricter)
# there is no quality filter, or possibiltiy to correct reads: plus for fastp
# --max-expected-errors could be used for quality filtering, expected errors from Phred Scores
# --max-n 0.1: filter out redas with more than 10% of N's



# ################################################################################
# # RUN FASTQC AFTER FASTP                                                       #
# ################################################################################



for FQ in "${OUTDIR_CUTADAPT}"/*.fastq.gz
do
    fastqc \
    -o ${OUTDIR_CUTADAPT}/fastqc/ \
    -t 32\
    ${FQ}
done


