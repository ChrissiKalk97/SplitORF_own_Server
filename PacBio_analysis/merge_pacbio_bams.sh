#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

pacbioRawDataDir=$1
pacbioMergedBamDir=$2

for dir in "${pacbioRawDataDir}"/*/; do
    bam1=$(find "$dir" -mindepth 2 -maxdepth 2 -name "*.bam" | sort | head -n 1)
    bam2=$(find "$dir" -mindepth 2 -maxdepth 2 -name "*.bam" | sort | tail -n 1)

    # echo "Merging:"
    # echo "$bam1"
    # echo "$bam2"
    # echo $(basename $dir)

    experiment=$(basename $dir)


    pbmerge -o ${pacbioMergedBamDir}/${experiment}_merged.bam $bam1 $bam2
done
  