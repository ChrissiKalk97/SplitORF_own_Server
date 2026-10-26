#!/bin/bash

# ----- This script analyses the HUVEC and CM Masspec Data from June 2025 ----- #
# ----- The PD results are first analyzed assumint that they are correct  ----- #
# ----- but actually the marking did not work, the second script performs ----- #
# ----- the tryptic digest of 2 possible references and really finds the  ----- #
# ----- unique peptides  ----- #


out_dir="/projects/splitorfs/work/Masspec/NewFasta_newAnalysis_23_10_2026/2026_10_23_analysis_results_with_ref"

data_dir="/projects/splitorfs/work/Masspec/NewFasta_newAnalysis_23_10_2026"

map_dir_old="/projects/splitorfs/work/split-orf-prediction/Output/tama_merged_iso_strin_mando_NMD_RI_masspec_files_20_07_26"

ms_input_dir="/projects/splitorfs/work/split-orf-prediction/Output/tama_merged_iso_strin_mando_NMD_RI_masspec_10000_10000_files_20_10_26"

so_path="/projects/splitorfs/work/split-orf-prediction/Output"


mkdir -p "$out_dir"
mkdir -p "$out_dir/background"
mkdir -p "$out_dir/coordinates"


eval "$(conda shell.bash hook)"
conda activate SplitORF

################################################################################
# MS ANALYSIS
################################################################################
# python peptide_groups_pd_analysis_with_reference_19_09_25.py \
#  --peptides_file "${data_dir}/20250430_AS_LC4_MAA_20049_01_VLD_HUVEC_F_PeptideGroups.txt" \
#  --so_id_mapping_file "${ms_input_dir}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_10000_10000_unique_SO_ID_files_20_10_26_so_id_mapping_with_assembly_info.tsv" \
#  --ref_id_mapping "${ms_input_dir}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_with_Uniprot_408_ref_unique_SO_ID_with_assembly_info_20_10_26_ref_id_mapping.tsv" \
#  --cell_type "huvec" \
#  --outdir "$out_dir" \
#  --uniprot_ensembl_mapping "${map_dir_old}/uniprot_ensembl_mapping.tsv"


# python peptide_groups_pd_analysis_with_reference_19_09_25.py \
#  --peptides_file "${data_dir}/20250505_AS_LC4_MAA_20050_01_VLD_iPSC_F_PeptideGroups.txt" \
#  --so_id_mapping_file "${ms_input_dir}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_10000_10000_unique_SO_ID_files_20_10_26_so_id_mapping_with_assembly_info.tsv" \
#  --ref_id_mapping "${ms_input_dir}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_with_Uniprot_408_ref_unique_SO_ID_with_assembly_info_20_10_26_ref_id_mapping.tsv" \
#  --cell_type "cm" \
#  --outdir "$out_dir" \
#  --uniprot_ensembl_mapping "${map_dir_old}/uniprot_ensembl_mapping.tsv"

################################################################################
# PREPARE UNIQUE PEPTIDES FOR RIBO-COV RUN (GENOMIC COORDS)
################################################################################
# for assembly in "NMD" "RI" "TAMA_HUVEC" "TAMA_CM"; do
#     python prepare_unique_peptides_for_genomic_positions.py \
#     --unique_peptide_information_csv "$out_dir"/huvec_validated_SO_protein_original_Ids_with_assembly.csv \
#     --favorite_assembly "$assembly" \
#     --cell_type HUVEC
# done

# for assembly in "NMD" "RI" "TAMA_HUVEC" "TAMA_CM"; do
#     python prepare_unique_peptides_for_genomic_positions.py \
#     --unique_peptide_information_csv "$out_dir"/cm_validated_SO_protein_original_Ids_with_assembly.csv \
#     --favorite_assembly "$assembly" \
#     --cell_type CM
# done


# Get all .sh files into an array
unique_pep_coord_files_huvec=(
    "${out_dir}/coordinates/TAMA_HUVEC_unique_peptides_of_HUVEC.bed" 
    "${out_dir}/coordinates/TAMA_CM_unique_peptides_of_HUVEC.bed" 
    "${out_dir}/coordinates/NMD_unique_peptides_of_HUVEC.bed" 
    "${out_dir}/coordinates/RI_unique_peptides_of_HUVEC.bed" 
    )

unique_pep_coord_files_cm=(
    "${out_dir}/coordinates/TAMA_HUVEC_unique_peptides_of_CM.bed" 
    "${out_dir}/coordinates/TAMA_CM_unique_peptides_of_CM.bed" 
    "${out_dir}/coordinates/NMD_unique_peptides_of_CM.bed" 
    "${out_dir}/coordinates/RI_unique_peptides_of_CM.bed" 
    )

coord_files=(
    "${so_path}/run_07.10.2026-11.40.34_HUVEC_ref_ORFanage_FIRST_and_Ens110/HUVECCoordsOfTranscriptsForSOConcatEns110.txt_transcript_positions.bed"
    "${so_path}/run_07.10.2026-10.43.25_CM_assembly_ref_ORFanage_FIRST_and_Ens110/CMCoordsOfTranscriptsForSOConcatEns110.txt_transcript_positions.bed"
    "${so_path}/run_07.04.2026-16.05.28_NMD_cont_subtraction/ExonCoordsWIthChr110_transcript_positions.bed"
    "${so_path}/run_07.04.2026-16.10.51_RI_contamination_subtraction/ExonCoordsWIthChr110_transcript_positions.bed"
    )



# get genomic coords huvec unique peptides
for i in "${!unique_pep_coord_files_huvec[@]}"; do

    unique_pep_file=${unique_pep_coord_files_huvec[$i]}
    coord_file=${coord_files[$i]}
    out_dir_huvec=$(dirname $unique_pep_file)
    outname=$(basename $unique_pep_file .bed)
    mkdir -p $out_dir_huvec/${outname}
    python /home/ckalk/scripts/SplitOrfs/split-orf-prediction/Genomic_scripts_18_10_24/genomic_dna_regions_polars.py \
     $unique_pep_file \
     $coord_file\
     "$out_dir_huvec/${outname}/Unique_DNA_Regions_genomic_final.bed"
done


# get genomic coords cm unique peptides
for i in "${!unique_pep_coord_files_cm[@]}"; do
    unique_pep_file=${unique_pep_coord_files_cm[$i]}
    coord_file=${coord_files[$i]}
    out_dir_cm=$(dirname $unique_pep_file)
    outname=$(basename $unique_pep_file .bed)
    mkdir -p $out_dir_cm/${outname}
    python /home/ckalk/scripts/SplitOrfs/split-orf-prediction/Genomic_scripts_18_10_24/genomic_dna_regions_polars.py \
     $unique_pep_file \
     $coord_file\
     "$out_dir_cm/${outname}/Unique_DNA_Regions_genomic_final.bed"
done

# concatenate all URs
mkdir -p "${out_dir}/coordinates/HUVEC_unique_peptides_concat"
mkdir -p "${out_dir}/coordinates/CM_unique_peptides_concat"

shopt -s globstar nullglob

out_cm="${out_dir}/coordinates/CM_unique_peptides_concat/Unique_DNA_Regions_genomic_final.bed"
> "$out_cm"                      # create or empty the output file for the concat URs
for bed_file in ${out_dir}/coordinates/*CM/*genomic_final.bed; do
    echo "Adding $bed_file"
    cat "$bed_file" >> "$out_cm"
done

out_huvec="${out_dir}/coordinates/HUVEC_unique_peptides_concat/Unique_DNA_Regions_genomic_final.bed"
> "$out_huvec"                      # create or empty the output file for the concat URs
for bed_file in ${out_dir}/coordinates/*HUVEC/*genomic_final.bed; do
    echo "Adding $bed_file"
    cat "$bed_file" >> "$out_huvec"
done

conda activate split-orf-pipeline-oct-26
ribo-cov /home/ckalk/scripts/SplitORFs/Masspec/Masspec_HUVEC_CM_2025/2026_10_26_ribo_cov_CM_IPSC_public_control_peptides.json

# after the conversion
# python helper_scripts/get_ensembl_gene_ids_from_uniprot_mapping.py \
#  --uniprot_ensembl_tsv ${OUTDIR}/background/Uniprot_to_ensembl_all_genes_HUVEC.tsv

# After GO analysis with SO genes (server)
# python /Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/gprofiler_scripts/create_grprofiler_barplot.py \
#  ${OUTDIR}/validated_genes_for_go/gProfiler_hsapiens_01-10-2025_15-26-41__intersections_high_set_expressed_proteins_genes_background_for_plotting.csv \
#  ${OUTDIR}/validated_genes_for_go/gprofiler_nice_barplot.svg \
# ";"