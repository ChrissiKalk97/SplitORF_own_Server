#!/bin/bash
. ~/spack/share/spack/setup-env.sh
spack load kallisto
read_list=("/scratch/fuchs/agschulz/kalk/sra_results/files/trimmed/SRR4081237_trimmomatic_crop_75_slow_2_threads.fastq"\
 "/scratch/fuchs/agschulz/kalk/sra_results/files/trimmed/SRR4081238_trimmomatic_crop_75_slow_2_threads.fastq"\
  "/scratch/fuchs/agschulz/kalk/sra_results/files/trimmed/SRR4081239_trimmomatic_crop_75_slow_2_threads.fastq"\
   "/scratch/fuchs/agschulz/kalk/sra_results/files/trimmed/SRR4081246_trimmomatic_crop_75_slow_2_threads.fastq"\
    "/scratch/fuchs/agschulz/kalk/sra_results/files/trimmed/SRR4081247_trimmomatic_crop_75_slow_2_threads.fastq"\
     "/scratch/fuchs/agschulz/kalk/sra_results/files/trimmed/SRR4081248_trimmomatic_crop_75_slow_2_threads.fastq")


for reads in "${read_list[@]}"; do
kallisto quant -i /scratch/fuchs/agschulz/kalk/Squanti3/kallisto/dKD_control/final/index/ref/reference\
 -o /scratch/fuchs/agschulz/kalk/Squanti3/kallisto/dKD_control/final/quant/ref/$(basename $reads _trimmomatic_crop_75_slow_2_threads.fastq)\
 -b 100\
 --single\
 -s 1\
 -l 75\
 $reads
done