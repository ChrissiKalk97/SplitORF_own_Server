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
# source preprocessing/fastqc_multiqc.sh ${INDIR} ${OUTDIR_FASTQC1} fastqc_unprocessed_multiqc


################################################################################
# RUN CUTADAPT TO TRIM ADAPTERS                                                #
################################################################################
OUTDIR_CUTADAPT=$3

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
# # RUN FASTQC AFTER CUTADAPT                                                    #
# ################################################################################
# source preprocessing/fastqc_multiqc.sh ${OUTDIR_CUTADAPT} ${OUTDIR_CUTADAPT} cutadapt_multiqc


################################################################################
# EXTRACT UMIS CUSTOM                                                          #
################################################################################
Umi_adpt_trimmed_path=$4
# python preprocessing/extract_umi/extract_compare_umis.py \
#  ${OUTDIR_CUTADAPT} \
#  ${Umi_adpt_trimmed_path}
# this is not optimal, think I prefer to to have the output in the same .out as the rest
# > preprocessing/extract_umi/extract_compare_umis_cutadapt.out 2>&1



# source preprocessing/fastqc_multiqc.sh ${Umi_adpt_trimmed_path} ${Umi_adpt_trimmed_path} UMI_trimmed_multiqc

cd ${Umi_adpt_trimmed_path}
gzip *.fastq
rm *.fastq

cd -

################################################################################
# FILTERING ETC FASTP                                                          #
################################################################################
fastpOut=$5
fastpFASTQC=$6
source preprocessing/filter_fastp.sh ${Umi_adpt_trimmed_path} ${fastpOut} cutadapt_umi_fastp
source preprocessing/fastqc_multiqc.sh ${fastpOut} ${fastpFASTQC} fastqc_fastp_trim_after_umi_trim_multiqc
multiqc --force --filename ${fastpOut}/fastp_filter_after_umi_trim_multiqc ${fastpOut}




