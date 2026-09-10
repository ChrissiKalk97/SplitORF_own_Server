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
conda activate pygtftk


for mode in "FIRST" "BEST"; do
    for cell_type in "CM" "HUVEC"; do
        if [[ "$cell_type" == "CM" ]]; then
            OUT_DIR=$OUT_DIR_CM
            mkdir -p $OUT_DIR_CM
            GTF=$CM_GTF
        else
            OUT_DIR=$OUT_DIR_HUVEC
            GTF=$HUVEC_GTF
        fi
        if [[ ! -d "$OUT_DIR/Orfanage_${mode}_09_09_26" ]]; then
            mkdir $OUT_DIR/Orfanage_${mode}_09_09_26
        fi

        # CDSs of genes that are present within the custom assembly
        python ${SCRIPT_DIR}/get_gtf_reference_prot_coding_trans.py \
        ${ENSEMBL_FILTERED_GTF} \
        ${GTF} \
        ${OUT_DIR}/Ens_110_prot_coding_filtered_CDS_for_${cell_type}_TAMA_09_09_26.gtf

        cd ${WORK_DIR}
        # the reference are all protein transcripts in Ensembl, not filtered
        ./orfanage --mode ${mode}  --reference ${GENOME_FASTA} \
        --query ${GTF}  \
        --output $OUT_DIR/Orfanage_${mode}_09_09_26/${cell_type}_TAMA_ORFanage_${mode}.gtf  ${OUT_DIR}/Ens_110_prot_coding_filtered_CDS_for_${cell_type}_TAMA_09_09_26.gtf


        python ${SCRIPT_DIR}/filter_CDS_entries.py \
        $OUT_DIR/Orfanage_${mode}_09_09_26/${cell_type}_TAMA_ORFanage_${mode}.gtf \
        $OUT_DIR/Orfanage_${mode}_09_09_26/${cell_type}_TAMA_ORFanage_${mode}_CDS.gtf

        python ${SCRIPT_DIR}/number_exons.py $OUT_DIR/Orfanage_${mode}_09_09_26/${cell_type}_TAMA_ORFanage_${mode}_CDS.gtf \
        $OUT_DIR/Orfanage_${mode}_09_09_26/${cell_type}_TAMA_ORFanage_${mode}_CDS_numbered.gtf
        

        bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_fiftynt_on_assembly.sh \
            $OUT_DIR/Orfanage_${mode}_09_09_26/${cell_type}_TAMA_ORFanage_${mode}_CDS_numbered.gtf \
            /home/ckalk/tools/NMD_fetaure_composition \
            $GENOME_FASTA \
            $ENSEMBL_FILTERED_GTF \
            ${cell_type}_TAMA_ORFanage_${mode}_09_09_26.csv
    done
done