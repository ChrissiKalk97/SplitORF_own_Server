#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq



################################################################################
# GET CONTAMINATION SEQS                                                       #
################################################################################
source preprocessing/contamination_filtering/prepare_Ignolia_cont_refs.sh
source preprocessing/construct_reference/concat_ignolia_reference.sh


fastpOut="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim"
Bowtie2_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
Bowtie2_base_name="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/index"
Bowtie2_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"

source alignments/bowtie2_align.sh ${Bowtie2_base_name} ${Bowtie2_ref_fasta} ${fastpOut} ${Bowtie2_out_dir} concat_transcriptome

