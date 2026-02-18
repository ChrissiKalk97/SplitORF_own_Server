#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Riboseq

# orfquant_path=$1
# ur_nmd=$2
# ur_ri=$3

orfquant_path="/projects/splitorfs/work/Riboseq/Output/ORFquant/ORFquant_output"
ur_nmd="/projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_NMD_16_12_24.bed"
ur_ri="/projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_RI_16_12_24.bed"


for gtf in ${orfquant_path}/*.gtf; do
    sample=${orfquant_path}/$(basename $gtf .gtf)
    gtf2bed < $gtf > ${sample}.bed
    bedtools intersect -wo -f 0.3 -a $ur_nmd -b ${sample}.bed > ${sample}_intersect_nmd.bed
    bedtools intersect -wo -f 0.3 -a $ur_ri -b ${sample}.bed > ${sample}_intersect_ri.bed
done

for bed in ${orfquant_path}/*_intersect_*.bed; do
    python intersection_ur_with_orfquant_orfs.py $bed
done