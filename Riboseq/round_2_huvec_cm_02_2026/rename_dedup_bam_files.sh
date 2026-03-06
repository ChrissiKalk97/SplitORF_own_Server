#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

dedup_bam_dir=$1


shopt -s nullglob

files=("${dedup_bam_dir}"/*.bam*)


for bam in "${files[@]}"
do
    filename=$(basename "$bam")
    # are the bamfiles genome or transcriptome aligned
    # remove respective suffix to get the sample name
    if [[ "$filename" == *".cutadapt_umi_fastp.bowtie1_concat_transcriptome_k1_R1_sorted_filtered_q10.bam"* ]]; then
        sample="${filename/.cutadapt_umi_fastp.bowtie1_concat_transcriptome_k1_R1_sorted_filtered_q10.bam/}"
        echo "${dedup_bam_dir}/${sample}"
        mv "${bam}" "${dedup_bam_dir}/${sample}"
    fi 
done

