#!/bin/bash

#----- This script performs the preprocessing steps for the Riboseq samples ----- #
# ----- first cutadapt is used to trim adapters and filter for quality ----- #
# ----- UMIs are extracted with a custom script as these are present on both reads  ----- #
# ----- Lastly, fastp is applied to correct bases and filter for a length range of 25-50 bp  ----- #
# ----- After each step FASTQC and MULTIQC reports are made  ----- #



eval "$(conda shell.bash hook)"
conda activate Riboseq

################################################################################
# RUN FASTQ                                                                    #
################################################################################
INDIR=$1
OUTDIR_FASTQC1=$2

# ################################################################################
# # RUN FASTQC                                                                   #
# ################################################################################
if [ ! -d ${OUTDIR_FASTQC1} ]; then
        mkdir ${OUTDIR_FASTQC1}
fi

source /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/preprocessing/fastqc_multiqc.sh \
 ${INDIR} ${OUTDIR_FASTQC1} fastqc_unprocessed_multiqc


################################################################################
# RUN CUTADAPT TO TRIM ADAPTERS                                                #
################################################################################
OUTDIR_CUTADAPT=$3

for FQ in "${INDIR}"/*R1.fastq.gz
do
SAMPLE=$(basename "$FQ" .R1.fastq.gz)
 cutadapt \
    --cores 0 \
    --pair-adapters \
    --quality-cutoff 20 \
    --adapter TGGAATTCTCGGGTGCCAAGG \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
    --output ${OUTDIR_CUTADAPT}/${SAMPLE}.cut.R1.fastq.gz \
    --paired-output ${OUTDIR_CUTADAPT}/${SAMPLE}.cut.R2.fastq.gz \
    --minimum-length 16 \
    ${FQ} \
    ${INDIR}/${SAMPLE}.R2.fastq.gz
done

# Parameters: as specified in the NEXTFLEX Small RNA-Seq Kit v4 instructions




# ################################################################################
# # RUN FASTQC AFTER CUTADAPT                                                    #
# ################################################################################
if [ ! -d ${OUTDIR_CUTADAPT} ]; then
        mkdir ${OUTDIR_CUTADAPT}  
fi


if [ ! -d ${OUTDIR_CUTADAPT}/fastqc ]; then
        mkdir ${OUTDIR_CUTADAPT}/fastqc
fi

source /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/preprocessing/fastqc_multiqc.sh \
 ${OUTDIR_CUTADAPT} ${OUTDIR_CUTADAPT}/fastqc cutadapt_multiqc





################################################################################
# FILTERING ETC FASTP                                                          #
################################################################################
fastpOut=$4
fastpFASTQC=$5
if [ ! -d ${fastpOut} ]; then
        mkdir ${fastpOut}
fi

if [ ! -d ${fastpFASTQC}  ]; then
        mkdir ${fastpFASTQC} 
fi

source /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/preprocessing/filter_fastp.sh \
 ${OUTDIR_CUTADAPT} ${fastpOut} fastp
source /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/preprocessing/fastqc_multiqc.sh \
 ${fastpOut} ${fastpFASTQC} fastqc_fastp_multiqc
multiqc --force --filename ${fastpOut}/fastp_multiqc ${fastpOut}




