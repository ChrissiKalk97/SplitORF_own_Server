#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq


transcriptome_assembly=$1 
genome_fasta=$2
transcriptome_fasta=$3
index_path=$4

gffread $transcriptome_assembly -g $genome_fasta -w $transcriptome_fasta



kallisto index -i ${index_path}.idx \
 $transcriptome_fasta