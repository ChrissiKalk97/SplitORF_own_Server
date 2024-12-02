#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate sra_tools
prefetch\
 --option-file /home/ckalk/scripts/SplitORFs/Riboseq/obtain_data/leukemia/SRR_Acc_List_leukemia_human.txt\
 -O /projects/splitorfs/work/Riboseq/data/leukemia/sra
