#!/bin/bash


# =============================================================================
# Script Name: 
# Description: This script runs ORFanage on the HUVEC TAMA assembly
# Usage:       bash 
# Author:      Christina Kalk
# Date:        2025-09-29
# =============================================================================

WORK_DIR="/home/ckalk/tools/ORFanage"

ENSEMBL_FILTERED_GTF="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
GENOME_FASTA="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
HUVEC_GTF="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/HUVEC/HUVEC_LR_SR_support_filtered.gtf"
CM_GTF="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/CM/CM_LR_SR_support_filtered.gtf"
OUT_DIR_HUVEC=$(dirname $HUVEC_GTF)/Orfanage
OUT_DIR_CM=$(dirname $CM_GTF)/Orfanage

SCRIPT_DIR="/home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts"

mkdir -p ${OUT_DIR_HUVEC}

eval "$(conda shell.bash hook)"
conda activate test-splitorf 


#################################################################################
# ------------------ CREATE REFERENCES FOR SPLIT-ORF PREDICTION --------------- #
#################################################################################
# # need to filter this for NMD transcripts and as well as for "protein coding transcripts"
# ~/tools/SplitORF_pipeline/Input2023/HUVEC_CM_assemblies/${cell_type}_10000_10000_tama_merged_assembly_transcriptome_gID_tID.fa
# ~/tools/SplitORF_pipeline/Input2023/HUVEC_CM_assemblies/${cell_type}_10000_10000_merged_tama_ExonCoordsOfTranscriptsForSO.txt

# # ORFanage CSV
# ~/tools/NMD_fetaure_composition/Output/CM_TAMA_ORFanage_FIRST_09_09_26/CM_TAMA_ORFanage_FIRST_09_09_26.csv
# # 50nt CSV
# ~/tools/NMD_fetaure_composition/Output/CM_merged_tama_10000_10000_iso_mando_stringtie_50nt/CM_merged_tama_10000_10000_iso_mando_stringtie_50nt.csv

# or does it make more sense to use the FASTA with all transcripts (so unfiltered, before LR and SR filtering) and then filter accordingly 
# and also change the header name accordingly?
# /projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/kallisto/${cell_type}_tama_merged_assembly_transcriptome.fa


for cell_type in "HUVEC" "CM"; do
    # from these two CSV files: get the transcripts that are NMD in either
    # as well as the transcripts with ORF that are NMD in neither (but predicted with ORFanage)

    # ------------------ get NMD and protein-coding transcript IDs --------------- #
    python get_nmd_prot_coding_transcripts.py \
    ~/tools/NMD_fetaure_composition/Output/${cell_type}_TAMA_ORFanage_FIRST_09_09_26/${cell_type}_TAMA_ORFanage_FIRST_09_09_26.csv \
    ~/tools/NMD_fetaure_composition/Output/${cell_type}_merged_tama_10000_10000_iso_mando_stringtie_50nt/${cell_type}_merged_tama_10000_10000_iso_mando_stringtie_50nt.csv \
    ${cell_type} \
    ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly

    # ------------------ obtain exon coords of custom assembly and concatenate with Ensembl exon coords --------------- #
    conda activate pygtftk
    # get exon coordinates
    python /home/ckalk/scripts/SplitOrfs/split-orf-prediction/Genomic_scripts_18_10_24/get_exon_coords_from_gtf.py \
    /projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/${cell_type}/${cell_type}_LR_SR_support_filtered.gtf \
    ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_10000_10000_merged_tama_ExonCoordsOfTranscriptsForSO.txt
    # concatenate all exon coordinates needed
    cat ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_10000_10000_merged_tama_ExonCoordsOfTranscriptsForSO.txt \
    ~/tools/SplitORF_pipeline/Input2023/ExonCoordsWIthChr110.bed > ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}CoordsOfTranscriptsForSOConcatEns110.txt

    # ------------------ get genomic CDS coords custom assembly concatenate with Ensembl genomic CDS coords --------------- #
    # get the CDS genomic coordinates
    python /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/get_CDS_genomic_coords_from_gtf.py \
    /projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/${cell_type}/Orfanage/Orfanage_FIRST_09_09_26/${cell_type}_TAMA_ORFanage_FIRST_CDS_numbered.gtf \
    ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/CDS_genomic_coords_${cell_type}_ORFanage_FIRST.bed

    # concatenate with Ensembl CDS genomic coordinates
    cat ~/tools/SplitORF_pipeline/Input2023/CDS_110_filtered_with_contaminants.bed \
    ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/CDS_genomic_coords_${cell_type}_ORFanage_FIRST.bed \
    > ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/CDS_genomic_coords_${cell_type}_assembly_Ens110_merged.bed

    # ------------------ get NMD and protein coding transcript sequences for SO pipeline input --------------- #
    # use the complete unfiltered (not SR, LR) assembly and just pick the selected 
    # transcripts
    # filter NMD transcripts
    conda activate Riboseq
    seqkit grep -f ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline.txt \
     /projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/kallisto/${cell_type}_tama_merged_assembly_transcriptome.fa \
      -o ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline.fa
    
    seqkit grep -f ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_protein_coding_ORFanage_FIRST.txt \
     /projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/kallisto/${cell_type}_tama_merged_assembly_transcriptome.fa \
      -o ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_protein_coding_ORFanage_FIRST.fa

    python /home/ckalk/scripts/SplitOrfs/split-orf-prediction/Input_scripts/change_fasta_header_custom_isoforms.py \
    /projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/${cell_type}/${cell_type}_LR_SR_support_filtered.gtf \
    ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_protein_coding_ORFanage_FIRST.fa \
    ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_protein_coding_ORFanage_FIRST_gID_tID.fa

    python /home/ckalk/scripts/SplitOrfs/split-orf-prediction/Input_scripts/change_fasta_header_custom_isoforms.py \
    /projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/${cell_type}/${cell_type}_LR_SR_support_filtered.gtf \
    ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline.fa \
    ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID.fa


    # ------------------ get protein sequences to use as reference proteins in SO pipeline  --------------- #
    # get the protein coding sequences from ORFanage FIRST
    # filter GTF for CDS features
    conda activate Riboseq
    # # -y writes the protein sequence
    gffread /projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/${cell_type}/Orfanage/Orfanage_FIRST_09_09_26/${cell_type}_TAMA_ORFanage_FIRST_CDS_numbered.gtf\
     -g $GENOME_FASTA -y ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_protein_coding_ORFanage_FIRST_peptide_sequences.fa

     # concat the protein coding sequences Ens110 and custom assembly
     cat ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_protein_coding_ORFanage_FIRST_peptide_sequences.fa \
     ~/tools/SplitORF_pipeline/Input2023/TSL_eq_filtered_29_09_25/protein_coding_peptide_sequences_tsl_eq_filtered_29_09_25.fa \
     > ~/tools/SplitORF_pipeline/Input2023/${cell_type}_assembly/${cell_type}_protein_coding_ORFanage_FIRST_and_Ens110_merged_peptide_sequences.fa
done


split-orf-prediction /home/ckalk/scripts/SplitORFs/PacBio_analysis/merge_stringtie_mando_isoquant/split_orf_pipeline_input_CM_10000_10000.json
split-orf-prediction /home/ckalk/scripts/SplitORFs/PacBio_analysis/merge_stringtie_mando_isoquant/split_orf_pipeline_input_HUVEC_10000_10000.json


