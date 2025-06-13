#!/bin/bash

#----- This script merges the short RNA-seq samples that are split ----- #
#----- Some samples are split because add-sequencing was performed  ----- #

raw_data_dir=$1
merged_data_dir=$2

if [[ ! -d "$merged_data_dir" ]]; then
  mkdir $merged_data_dir
fi

for dir in ${raw_data_dir}/*/; do 
# this is not working need to find a way to combine corresponding files of each subdir into the 
# merged dir also write single files into the merged?

  fastq_array=($dir/*.fq.gz)
  nr_files=${#fastq_array[@]}
  sample_name=$(basename "${fastq_array[0]}")
  sample_name="${sample_name%_*_*_*_*}"

  R_1_file=${merged_data_dir}/${sample_name}_merged_1.fq.gz
  R_2_file=${merged_data_dir}/${sample_name}_merged_2.fq.gz

  if [[ ! -f "$R_1_file" ]]; then
    touch "$R_1_file"
    touch "$R_2_file"
  fi

  for file in "${fastq_array[@]}"; do
    if [[ $file == *_1.fq.gz ]]; then
      cat "$file" >> "$R_1_file"
    else
      cat "$file" >> "$R_2_file"
    fi

  done
  
  
done
