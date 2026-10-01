#!/bin/bash


# =============================================================================
# Script Name: 
# Description: This script runs ORFanage on the HUVEC TAMA assembly
# Usage:       bash 
# Author:      Christina Kalk
# Date:        2025-09-29
# =============================================================================

WORK_DIR="/home/ckalk/tools/ORFanage"

ENSEMBL_FULL_GTF="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"
GENOME_FASTA="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
HUVEC_GTF="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama/HUVEC/HUVEC_merged_tama_gene_id.gtf"
OUT_DIR=$(dirname $HUVEC_GTF)/Orfanage

SCRIPT_DIR="/home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts"


if [[ ! -d "$OUT_DIR" ]]; then
    mkdir $OUT_DIR
fi


eval "$(conda shell.bash hook)"
conda activate pygtftk
python ${SCRIPT_DIR}/get_gtf_reference_prot_coding_trans.py \
 ${ENSEMBL_FULL_GTF} \
 ${HUVEC_GTF} \
 ${OUT_DIR}/Ens_110_prot_coding_CDS_for_HUVEC_TAMA_29_09_25.gtf


if [[ ! -d "$OUT_DIR/HUVEC_TAMAv1_FIRST_29_09_25" ]]; then
    mkdir $OUT_DIR/HUVEC_TAMAv1_FIRST_29_09_25
fi

cd ${WORK_DIR}
# the reference are all protein transcripts in Ensembl, not filtered
./orfanage --mode FIRST  --reference ${GENOME_FASTA} \
--query ${HUVEC_GTF}  \
--output ${OUT_DIR}/HUVEC_TAMAv1_FIRST_29_09_25/HUVEC_TAMA_ORFanage_FIRST.gtf  ${OUT_DIR}/Ens_110_prot_coding_CDS_for_HUVEC_TAMA_29_09_25.gtf


python ${SCRIPT_DIR}/filter_CDS_entries.py \
 ${OUT_DIR}/HUVEC_TAMAv1_FIRST_29_09_25/HUVEC_TAMA_ORFanage_FIRST.gtf \
  ${OUT_DIR}/HUVEC_TAMAv1_FIRST_29_09_25/HUVEC_TAMA_ORFanage_FIRST_CDS.gtf

 python ${SCRIPT_DIR}/number_exons.py ${OUT_DIR}/HUVEC_TAMAv1_FIRST_29_09_25/HUVEC_TAMA_ORFanage_FIRST_CDS.gtf \
  ${OUT_DIR}/HUVEC_TAMAv1_FIRST_29_09_25/HUVEC_TAMA_ORFanage_FIRST_CDS_numbered.gtf
 

 bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_fiftynt_on_assembly.sh \
    ${OUT_DIR}/HUVEC_TAMAv1_FIRST_29_09_25/HUVEC_TAMA_ORFanage_FIRST_CDS_numbered.gtf \
    /home/ckalk/tools/NMD_fetaure_composition \
    $GENOME_FASTA \
    $ENSEMBL_FULL_GTF \
    TAMA_HUVEC_ORFanage_v1_29_09_25.csv