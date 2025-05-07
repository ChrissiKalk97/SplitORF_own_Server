#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq



################################################################################
# GET CONTAMINATION SEQS                                                       #
################################################################################
# source preprocessing/contamination_filtering/prepare_Ignolia_cont_refs.sh
# source preprocessing/construct_reference/concat_ignolia_reference.sh


fastpOut=$1
Bowtie2_ref_fasta=$2
Bowtie2_base_name=$3
Bowtie2_out_dir=$4

echo $fastpOut 
echo $Bowtie2_ref_fasta
echo $Bowtie2_base_name
echo $Bowtie2_out_dir

source alignments/bowtie2_align_k1.sh ${Bowtie2_base_name} no_index ${fastpOut} ${Bowtie2_out_dir} concat_transcriptome



