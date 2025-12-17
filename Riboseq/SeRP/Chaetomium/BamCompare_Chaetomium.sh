#!/bin/bash

bam_dir=$1
input_or_e_denominator=$2
input_or_e_numerator=$3
w_or_s_denominator=$4
w_or_s_numerator=$5
outdir=$6


eval "$(conda shell.bash hook)"
conda activate Riboseq
export TMPDIR=/scratch/tmp/$USER
for bam in $bam_dir/*.bam; do
  numerator_filename=$(basename "$bam")
  
  if [[ "$numerator_filename" =~ ^uf_muellermcnicoll_([0-9_]+)_MD_${input_or_e_numerator}_${w_or_s_numerator}([0-9]+)\.cut\.fastp\.bowtie1_Chaetomium_transcriptome_k1_R1_sorted_filtered_q10\.bam$ ]]; then
        date="${BASH_REMATCH[1]}"
        batch="${BASH_REMATCH[2]}"

                        
        denominator_filename=$bam_dir/uf_muellermcnicoll_2025_05_??_MD_${input_or_e_denominator}_${w_or_s_denominator}${batch}.cut.fastp.bowtie1_Chaetomium_transcriptome_k1_R1_sorted_filtered_q10.bam
        echo $numerator_filename
        echo $denominator_filename
        
        bamCompare -b1 $bam_dir/$numerator_filename\
         -b2 $denominator_filename \
         -o $outdir/${w_or_s_numerator}_${input_or_e_numerator}_over_${w_or_s_denominator}_${input_or_e_denominator}_${batch}_whole_trans_b1_no_smooth.bw \
         --operation ratio \
         --binSize 1 \
         -of bigwig\
         -p 64

  fi
done



