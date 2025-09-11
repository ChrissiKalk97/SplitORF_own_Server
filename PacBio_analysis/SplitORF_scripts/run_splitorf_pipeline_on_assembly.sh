#!/bin/bash

#----- This script performs the preprocessing steps and runs the SplitORF pipeline ----- #
# ----- for a custom assembly       ----- #

custom_gtf=$1
splitorf_pipeline=$2
genome_fasta=$3

custom_gtf_dir=$(dirname $custom_gtf)
custom_gtf_base=$(basename $custom_gtf .gtf)

eval "$(conda shell.bash hook)"
conda activate pacbio


#################################################################################
# ------------------ PREPROCESSING FOR SO PIPELINE           ------------------ #
#################################################################################
transcriptome_fasta=$custom_gtf_dir/${custom_gtf_base}.fasta
transcriptome_fasta_with_gid=$custom_gtf_dir/$(basename $transcriptome_fasta .fasta )withGid.fasta

# get transcriptome FASTA
# gffread $custom_gtf -g $genome_fasta -w $transcriptome_fasta

# # change the FASTA header to gene_id|transcript_id
# conda activate pygtftk

# python \
#  $splitorf_pipeline/Input_scripts/change_fasta_header_custom_isoforms.py \
#  $custom_gtf \
#  $transcriptome_fasta \
#  $transcriptome_fasta_with_gid


# python $splitorf_pipeline/Genomic_scripts_18_10_24/get_exon_coords_from_gtf.py \
#  $custom_gtf \
#  $custom_gtf_dir/${custom_gtf_base}_ExonCoordsOfTranscriptsForSO.txt


conda activate SplitORF
cd $splitorf_pipeline
bash $splitorf_pipeline/run_splitorfs_pipeline.sh \
 $splitorf_pipeline/Input2023/TSL_filtered/protein_coding_peptide_sequences_110_tsl_refseq_filtered.fa \
 $transcriptome_fasta_with_gid \
 $splitorf_pipeline//Input2023/ENSEMBLhuman110PFAMformat.bed \
 $splitorf_pipeline/Input2023/TSL_filtered/protein_coding_transcript_and_gene_cDNA_110_tsl_refseq_filtered.fa \
 $custom_gtf_dir/${custom_gtf_base}_ExonCoordsOfTranscriptsForSO.txt \
 blast