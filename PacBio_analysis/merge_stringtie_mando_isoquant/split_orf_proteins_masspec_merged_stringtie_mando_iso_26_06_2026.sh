#!/bin/bash

#----- This script runs gffcompare on 2 assemblies to merge and compare them ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio

so_path="/projects/splitorfs/work/split-orf-prediction/Output"
script_path="/home/ckalk/scripts/SplitOrfs/split-orf-prediction/Output_scripts/compare_ORFs_amongst_transcripts"

tama_so_cm_dir="run_25.06.2026-11.22.12_CM_mando_iso_stringtie_rescued"
tama_so_huvec_dir="run_25.06.2026-13.59.59_HUVEC_mando_iso_stringtie_rescued"
RI_so_dir="run_07.04.2026-16.05.28_NMD_cont_subtraction"
NMD_so_dir="run_07.04.2026-16.10.51_RI_contamination_subtraction"
masspec_dir="tama_merged_iso_strin_mando_NMD_RI_masspec_files_26_06_26"

if [[ ! -d "${so_path}/${masspec_dir}" ]]; then
    mkdir ${so_path}/${masspec_dir}
fi

python ${script_path}/get_unique_ids_for_orfs_for_masspec_assembly_names.py \
 --fasta3 "${so_path}/${tama_so_cm_dir}/Proteins_for_masspec.fa" \
 --fasta4 "${so_path}/${tama_so_huvec_dir}/Proteins_for_masspec.fa" \
 --assembly_list "NMD,RI,TAMA_CM,TAMA_HUVEC" \
 "${so_path}/${NMD_so_dir}/Proteins_for_masspec.fa" \
 "${so_path}/${RI_so_dir}/Proteins_for_masspec.fa" \
 "${so_path}/${masspec_dir}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_unique_SO_ID_26_06_26.fasta" \
 "SplitOrfProtein"


python ${script_path}/merge_reference_splitorf_prots_for_masspec.py \
 "/projects/splitorfs/work/reference_files/Homo_sapiens_sp_incl_isoforms_TaxID_9606_Release_408.fasta" \
 "${so_path}/${masspec_dir}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_unique_SO_ID_26_06_26.fasta" \
 "${so_path}/${masspec_dir}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_with_Uniprot_408_ref_unique_SO_ID_with_assembly_info_26_06_26.fasta"


 