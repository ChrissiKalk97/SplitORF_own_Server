#!/bin/bash

data_dir=$1
fastqc_dir=$2
multiQC_outname=$3
raw=$4

if [ ! -d ${fastqc_dir} ]; then
    mkdir $fastqc_dir
fi

if [ "$raw" == "raw" ];
then
    # change into the raw data directory
    cd $data_dir

    # empty array for FQ files
    fq_files=()
    while IFS= read -r line; do
        fq_files+=("$line")
    done < <(find "$(pwd)" -type f -name "*.fq.gz")
    # < (some command) is a process substitution in bash
    # execute command, treat output like a file


    cd -

    # how many files do we have?
    echo "${#fq_files[@]}"

else
    shopt -s nullglob
    fq_files=("${data_dir}"/*fastq.gz)
fi

for FQ in "${fq_files[@]}"; 
do
    (
        fastqc \
        -o ${fastqc_dir}/ \
        -t 32\
        ${FQ}
    )&
done

wait

multiqc --force --filename ${fastqc_dir}/${multiQC_outname} ${fastqc_dir}