#!/bin/bash

bam_dir=$1
chx_or_puro=$2
input_or_mock=$3
E=$4


eval "$(conda shell.bash hook)"
conda activate Riboseq
export TMPDIR=/scratch/tmp/$USER
for bam in $bam_dir/*.bam; do
  importin_filename=$(basename "$bam")
  if [[ "$importin_filename" =~ ^uf_muellermcnicoll_([0-9_]+)_RR_([AB][12])_${chx_or_puro}_E([0-9]+)\.cut\.fastp\.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10\.bam$ ]]; then
        date="${BASH_REMATCH[1]}"
        importin="${BASH_REMATCH[2]}"
        batch="${BASH_REMATCH[3]}"


        input_filename=$bam_dir/uf_muellermcnicoll_2025_05_??_RR_${input_or_mock}_${chx_or_puro}_${E}${batch}.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10.bam
        echo $input_filename
        echo $importin_filename
        
        bamCompare -b1 $bam_dir/$importin_filename\
         -b2 $input_filename\
         -o $bam_dir/enrichment_plots_CDS/whole_transcript_bigwig/${importin}_E${batch}_over_${input_or_mock}_${chx_or_puro}_whole_trans_b1_no_smooth.bw\
         --operation ratio \
         --binSize 1 \
         -of bigwig\
         -p 64

  fi
done



