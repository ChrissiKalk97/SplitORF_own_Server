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
umi_adpt_trimmed_path=$4
fastp_out=$5
fastp_fastqc=$6
module_dir=$7

cd "${module_dir}"

# ################################################################################
# # RUN FASTQC                                                                   #
# ################################################################################
source preprocessing/fastqc_multiqc.sh ${indir} ${outdir_fastqc1} fastqc_unprocessed_multiqc


################################################################################
# RUN CUTADAPT TO TRIM ADAPTERS                                                #
################################################################################
for fq in "${indir}"/*R1.fastq.gz
do
sample=$(basename "$fq" .R1.fastq.gz)
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o ${outdir_cutadapt}/${sample}.R1.cutadapt.fastq.gz \
    -p ${outdir_cutadapt}/${sample}.R2.cutadapt.fastq.gz \
    --cores 16 \
    --minimum-length 20 \
    --max-n 0.1 \
    --max-expected-errors 1 \
    ${fq} \
    ${indir}/${sample}.R2.fastq.gz
done

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
source preprocessing/fastqc_multiqc.sh ${outdir_cutadapt} ${outdir_cutadapt} cutadapt_multiqc


################################################################################
# EXTRACT UMIS CUSTOM                                                          #
################################################################################
python preprocessing/extract_umi/extract_compare_umis.py \
 ${outdir_cutadapt} \
 ${umi_adpt_trimmed_path}
# this is not optimal, think I prefer to to have the output in the same .out as the rest
# > preprocessing/extract_umi/extract_compare_umis_cutadapt.out 2>&1



source preprocessing/fastqc_multiqc.sh ${umi_adpt_trimmed_path} ${umi_adpt_trimmed_path} UMI_trimmed_multiqc

cd ${umi_adpt_trimmed_path}
gzip *.fastq
rm *.fastq

cd -

################################################################################
# FILTERING ETC FASTP                                                          #
################################################################################
source preprocessing/filter_fastp.sh ${umi_adpt_trimmed_path} ${fastp_out} cutadapt_umi_fastp
source preprocessing/fastqc_multiqc.sh ${fastp_out} ${fastp_fastqc} fastqc_fastp_trim_after_umi_trim_multiqc
multiqc --force --filename ${fastp_out}/fastp_filter_after_umi_trim_multiqc ${fastp_out}




