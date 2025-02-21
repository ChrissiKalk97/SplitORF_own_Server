#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq

python prepare_URs_for_Ribotricer.py\
 /projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_NMD_16_12_24.bed\
 /projects/splitorfs/work/Riboseq/Output/Ribotricer/URs_as_ORFs/NMD_URs_as_ORFs.tsv