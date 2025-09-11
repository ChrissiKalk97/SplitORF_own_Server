#!/bin/bash

#----- This script performs the preprocessing steps and runs the SplitORF pipeline ----- #
# ----- for a custom assembly       ----- #

custom_gtf=$1
fiftynt_pipeline=$2
genome_fasta=$3
ensembl_full_gtf=$4
outname=$5

export TMPDIR=/scratch/tmp/$USER

custom_gtf_dir=$(dirname $custom_gtf)
custom_gtf_base=$(basename $custom_gtf .gtf)

eval "$(conda shell.bash hook)"
conda activate pygtftk

cd $fiftynt_pipeline
python fifty_nt_rule.py \
    $custom_gtf \
    $ensembl_full_gtf \
    $genome_fasta \
    $outname



