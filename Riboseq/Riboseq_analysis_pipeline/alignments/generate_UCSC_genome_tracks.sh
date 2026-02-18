#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Riboseq


bam_dir=$1
bam_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"

for bam in ${bam_dir}/*_dedup_filtered.bam
do
    sample=$(basename ${bam} .bam)
    out_dir=$(dirname ${bam})
    bamCoverage -b ${bam} -o ${out_dir}/${sample}.bw \
    --normalizeUsing CPM --binSize 10 -p 32
done