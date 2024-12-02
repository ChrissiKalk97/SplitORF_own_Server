#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate sra_tools
for i in {08..11}; do

fasterq-dump /projects/splitorfs/work/Riboseq/data/leukemia/sra/SRR112946$i\
 -O /projects/splitorfs/work/Riboseq/data/leukemia/fastq

done
