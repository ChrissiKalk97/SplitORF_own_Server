#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq

data_dir=$1

bam_files=(${data_dir}/*.bam)

if [ ! -d ${data_dir}/fastq ];then
	mkdir ${data_dir}/fastq
fi

if [ ! -d ${data_dir}/fastq/fastqc ];then
	mkdir ${data_dir}/fastq/fastqc
fi

for bam in "${bam_files[@]}"; do

    sample=$(basename $bam .bam)

    if [ ! -e "${data_dir}/${sample}_name_sorted.bam" ]; then
        samtools sort -n ${bam} -o ${data_dir}/${sample}_name_sorted.bam
    fi

    if [ ! -e "${data_dir}/fastq/${sample}.fastq" ]; then
        bedtools bamtofastq -i ${data_dir}/${sample}_name_sorted.bam -fq ${data_dir}/fastq/${sample}.fastq 
    fi

done

source /home/ckalk/scripts/SplitORFs/UPF1_deletion/fastqc_multiqc_for_all.sh \
${data_dir}/fastq ${data_dir}/fastq/fastqc fastqc_from_bam_multiqc FASTQ