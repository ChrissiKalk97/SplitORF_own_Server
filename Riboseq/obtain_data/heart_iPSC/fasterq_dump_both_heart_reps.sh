#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate sra_tools
fasterq-dump /projects/splitorfs/work/Riboseq/data/heart_iPSC/sra/ERR3367798\
 -O /projects/splitorfs/work/Riboseq/data/heart_iPSC/fastq

fasterq-dump /projects/splitorfs/work/Riboseq/data/heart_iPSC/sra/ERR3367797\
 -O /projects/splitorfs/work/Riboseq/data/heart_iPSC/fastq
