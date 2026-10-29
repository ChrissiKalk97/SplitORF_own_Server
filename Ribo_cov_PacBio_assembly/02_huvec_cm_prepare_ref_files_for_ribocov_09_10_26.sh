#!/bin/bash


# =============================================================================
# Script Name: 02_huvec_cm_prepare_ref_files_for_ribocov_09_10_26.sh
# Description: This script generates the Input files needed for the ribo-cov
#               run with a custom assembly:
#              - Step 1: 
#              - Step 2: 
#              - Step 3: 
#              - Step 4: 
#              - Step 5: 
# Usage:       bash 02_huvec_cm_prepare_ref_files_for_ribocov_09_10_26.sh
# Author:      Christina Wrobel
# Date:        2026-10-09
# =============================================================================






eval "$(conda shell.bash hook)"
conda activate Riboseq

################################################################################
# PATH DEFINTIONS                                                              #
################################################################################



input_dir="/home/ckalk/tools/SplitORF_pipeline/Input2023"
input_dir_cm="${input_dir}"/CM_assembly/ribocov_input
input_dir_huvec="${input_dir}"/HUVEC_assembly/ribocov_input


mkdir -p "${input_dir_cm}"
mkdir -p "${input_dir_huvec}"

out_dir="/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2_PacBio_assembly"
out_dir_cm="${out_dir}"/CM
out_dir_huvec="${out_dir}"/HUVEC

mkdir -p "$out_dir_cm"
mkdir -p "$out_dir_huvec"

report_dir="/home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/outreports_of_runs"
report_dir_cm="${report_dir}"/CM
report_dir_huvec="${report_dir}"/HUVEC

mkdir -p "$report_dir_cm"
mkdir -p "$report_dir_huvec"

script_dir_02="/home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/02_huvec_cm_prepare_ref_files_for_ribocov"

# 1-based
cm_trans_coord="/home/ckalk/tools/SplitORF_pipeline/Input2023/CM_assembly/CMCoordsOfTranscriptsForSOConcatEns110.txt"
# 0-based
cm_cds_coords="/home/ckalk/tools/SplitORF_pipeline/Input2023/CM_assembly/CDS_genomic_coords_CM_assembly_Ens110_merged.bed"
cm_ur_dir="/projects/splitorfs/work/split-orf-prediction/Output/run_07.10.2026-10.43.25_CM_assembly_ref_ORFanage_FIRST_and_Ens110"

# 1-based
huvec_trans_coord="/home/ckalk/tools/SplitORF_pipeline/Input2023/HUVEC_assembly/HUVECCoordsOfTranscriptsForSOConcatEns110.txt"
# 0-based
huvec_cds_coords="/home/ckalk/tools/SplitORF_pipeline/Input2023/HUVEC_assembly/CDS_genomic_coords_HUVEC_assembly_Ens110_merged.bed"
huvec_ur_dir="/projects/splitorfs/work/split-orf-prediction/Output/run_07.10.2026-11.40.34_HUVEC_ref_ORFanage_FIRST_and_Ens110"

cds_coords_ensembl="/home/ckalk/tools/SplitORF_pipeline/Input2023/TSL_eq_filtered_29_09_25/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered.bed"

cds_coords_cm="/home/ckalk/tools/SplitORF_pipeline/Input2023/CM_assembly/CDS_genomic_coords_CM_ORFanage_FIRST.bed"
cds_coords_huvec="/home/ckalk/tools/SplitORF_pipeline/Input2023/HUVEC_assembly/CDS_genomic_coords_HUVEC_ORFanage_FIRST.bed"


######### CREATE THE 3' UTR GENOMIC COORDS WITHOUT CDS OVERLAP #################
######### AND WITHOUT UNIQUE REGION OVERLAP                    #################

bash "$script_dir_02"/get_3prime_genomic_coords.sh \
    -t "$cm_trans_coord" \
    -c "$cm_cds_coords" \
    -s "$script_dir_02" \
    -r "CM" \
    -o "$input_dir_cm" \
    -u "$cm_ur_dir"


bash "$script_dir_02"/get_3prime_genomic_coords.sh \
    -t "$huvec_trans_coord" \
    -c "$huvec_cds_coords" \
    -s "$script_dir_02" \
    -r "HUVEC" \
    -o "$input_dir_huvec" \
    -u "$huvec_ur_dir"


cat "$cds_coords_huvec" "$cds_coords_ensembl" > "$input_dir_huvec"/CDS_genomic_coords_HUVEC_assembly_Ens110_no_contamination.bed
cat "$cds_coords_cm" "$cds_coords_ensembl" > "$input_dir_cm"/CDS_genomic_coords_CM_assembly_Ens110_no_contamination.bed
