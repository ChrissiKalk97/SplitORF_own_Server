#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate sra_tools
prefetch\
 --option-file /home/ckalk/scripts/SplitORFs/Riboseq/obtain_data/heart_iPSC/SRR_Acc_List_only_human_heart_2_reps.txt\
 -O /projects/splitorfs/work/Riboseq/data/heart_iPSC/sra
