#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

bam_dir=$1
lima_outdir=$2
primer_fasta=$3


if [[ ! -d "$lima_outdir" ]]; then
    mkdir $lima_outdir
fi

shopt -s nullglob
bam_files=("${bam_dir}"/*bam)

for bam in "${bam_files[@]}"; 
do
    sample=$(basename $bam .bam)
    lima $bam $primer_fasta $lima_outdir/${sample}_lima.bam --isoseq 
done