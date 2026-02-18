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
indir=$1
outdir_fastqc1=$2
outdir_cutadapt=$3
fastp_dir=$4
fastp_fastqc_dir=$5
module_dir=$6

# ################################################################################
# # RUN FASTQC                                                                   #
# ################################################################################
if [ ! -d ${outdir_fastqc1} ]; then
        mkdir ${outdir_fastqc1}
fi

source "${module_dir}"/preprocessing/fastqc_multiqc.sh \
 ${indir} ${outdir_fastqc1} fastqc_unprocessed_multiqc


################################################################################
# RUN CUTADAPT TO TRIM ADAPTERS                                                #
################################################################################


for FQ in "${indir}"/*R1.fastq.gz
do
SAMPLE=$(basename "$FQ" .R1.fastq.gz)
 cutadapt \
    --cores 0 \
    --pair-adapters \
    --quality-cutoff 20 \
    --adapter TGGAATTCTCGGGTGCCAAGG \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
    --output ${outdir_cutadapt}/${SAMPLE}.cut.R1.fastq.gz \
    --paired-output ${outdir_cutadapt}/${SAMPLE}.cut.R2.fastq.gz \
    --minimum-length 16 \
    ${FQ} \
    ${indir}/${SAMPLE}.R2.fastq.gz
done

# Parameters: as specified in the NEXTFLEX Small RNA-Seq Kit v4 instructions




# ################################################################################
# # RUN FASTQC AFTER CUTADAPT                                                    #
# ################################################################################
if [ ! -d ${outdir_cutadapt} ]; then
        mkdir ${outdir_cutadapt}  
fi


if [ ! -d ${outdir_cutadapt}/fastqc ]; then
        mkdir ${outdir_cutadapt}/fastqc
fi

source "${module_dir}"/preprocessing/fastqc_multiqc.sh \
 ${outdir_cutadapt} ${outdir_cutadapt}/fastqc cutadapt_multiqc





################################################################################
# FILTERING ETC FASTP                                                          #
################################################################################

if [ ! -d ${fastp_dir} ]; then
        mkdir ${fastp_dir}
fi

if [ ! -d ${fastp_fastqc_dir}  ]; then
        mkdir ${fastp_fastqc_dir} 
fi

source "${module_dir}"/preprocessing/filter_fastp.sh \
 ${outdir_cutadapt} ${fastp_dir} fastp
source "${module_dir}"/preprocessing/fastqc_multiqc.sh \
 ${fastp_dir} ${fastp_fastqc_dir} fastqc_fastp_multiqc
multiqc --force --filename ${fastp_dir}/fastp_multiqc ${fastp_dir}




