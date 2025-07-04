#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate longqc-env

bam_dir=$1
longqc_dir=$2

if [[ ! -d "$longqc_dir" ]]; then
    mkdir $longqc_dir
fi

if [[ ! -d "$longqc_dir"/sampleqc ]]; then
    mkdir $longqc_dir/sampleqc
fi

if [[ ! -d "$longqc_dir"/runqc ]]; then
    mkdir $longqc_dir/runqc
fi

shopt -s nullglob
bam_files=("${bam_dir}"/*bam)

cd /home/ckalk/scripts/SplitORFs/PacBio_analysis/LongQC
for bam in "${bam_files[@]}"; 
do
    ./longQC.py sampleqc -x pb-sequel -o $longqc_dir/sampleqc $bam
    ./longQC.py runqc -x pb-sequel -o $longqc_dir/runqc $bam
done