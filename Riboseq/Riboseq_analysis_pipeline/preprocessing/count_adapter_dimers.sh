#!/bin/bash

raw_file_directory="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_1"

for fastq in $raw_file_directory/uf*R1.fastq.gz
do
    sample=$(basename $fastq R1.fastq.gz)
    nr_full_adapter=$(zgrep -c "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" $fastq)
    echo $sample
    echo $nr_full_adapter
    nr_9N_beginning=$(zgrep -c "^NNNNNNNNN" $fastq)
    echo $nr_9N_beginning
done

python analyze_adapter_dimers.py $raw_file_directory