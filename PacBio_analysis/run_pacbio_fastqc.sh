#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

bam_dir=$1
fastqc_dir=$2
multiqc_outname=$3

if [[ ! -d "$fastqc_dir" ]]; then
    mkdir $fastqc_dir
fi


shopt -s nullglob
bam_files=("${bam_dir}"/*bam)


# for bam in "${bam_files[@]}"; 
# do
#     (
#         fastqc \
#         -o ${fastqc_dir}/ \
#         -t 32\
#         ${bam}
#     )&
# done

# wait

multiqc --force --filename ${fastqc_dir}/${multiqc_outname} ${fastqc_dir}