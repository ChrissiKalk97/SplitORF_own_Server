#!/bin/bash

#----- This script runs gffcompare on 2 assemblies to merge and compare them ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio

outdir_tama="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama"
so_path="/home/ckalk/tools/SplitORF_pipeline/Output"
script_path="${so_path}/compare_ORFs_amongst_transcripts"

tama_so_cm_dir="run_12.09.2025-14.10.14_CM_tama_merged"
tama_so_huvec_dir="run_12.09.2025-17.51.04_HUVEC_tama_merged"
RI_so_dir="run_18.06.2025-11.59.36_RI_transcripts"
NMD_so_dir="run_18.06.2025-09.35.29_NMD_transcripts"
masspec_dir="tama_NMD_RI_masspec_files"

if [[ ! -d "${so_path}/${masspec_dir}" ]]; then
    mkdir ${so_path}/${masspec_dir}
fi

python ${script_path}/get_unique_ids_for_orfs_for_masspec_assembly_names.py \
 --fasta3 "${so_path}/${tama_so_cm_dir}/Proteins_for_masspec.fa" \
 --fasta4 "${so_path}/${tama_so_huvec_dir}/Proteins_for_masspec.fa" \
 --assembly_list "NMD,RI,TAMA_CM,TAMA_HUVEC" \
 "${so_path}/${NMD_so_dir}/Proteins_for_masspec.fa" \
 "${so_path}/${RI_so_dir}/Proteins_for_masspec.fa" \
 "${so_path}/${masspec_dir}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_unique_SO_ID_16_09_25.fasta" \
 "SplitOrfProtein"


python ${script_path}/merge_reference_splitorf_prots_for_masspec.py \
 "/projects/splitorfs/work/reference_files/Homo_sapiens_sp_incl_isoforms_TaxID_9606_Release_408.fasta" \
 "/home/ckalk/tools/SplitORF_pipeline/Output/tama_NMD_RI_masspec_files/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_unique_SO_ID_16_09_25.fasta" \
 "/home/ckalk/tools/SplitORF_pipeline/Output/tama_NMD_RI_masspec_files/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_with_Uniprot_408_ref_unique_SO_ID_with_assembly_info.fasta"


 