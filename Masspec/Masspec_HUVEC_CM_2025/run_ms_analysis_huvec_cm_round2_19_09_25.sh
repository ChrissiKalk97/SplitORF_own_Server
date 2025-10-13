#!/bin/bash

# ----- This script analyses the HUVEC and CM Masspec Data from June 2025 ----- #
# ----- The PD results are first analyzed assumint that they are correct  ----- #
# ----- but actually the marking did not work, the second script performs ----- #
# ----- the tryptic digest of 2 possible references and really finds the  ----- #
# ----- unique peptides  ----- #


OUTDIR="/projects/splitorfs/work/Masspec/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25"

DATADIR="/projects/splitorfs/work/Masspec/New_MS_run_19_09_25_tama_assembly_SOs"

MAPDIR="/home/ckalk/tools/SplitORF_pipeline/Output/tama_NMD_RI_masspec_files"

SO_PATH="/home/ckalk/tools/SplitORF_pipeline/Output"

if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

if [ ! -d "$OUTDIR/background" ]; then
    mkdir -p "$OUTDIR/background"
fi

if [ ! -d "$OUTDIR/coordinates" ]; then
    mkdir -p "$OUTDIR/coordinates"
fi

eval "$(conda shell.bash hook)"
conda activate SplitORF

python peptide_groups_pd_analysis_with_reference_19_09_25.py \
 --peptides_file "${DATADIR}/20250430_AS_LC4_MAA_20049_01_VLD_HUVEC_F_PeptideGroups.txt" \
 --so_id_mapping_file "${MAPDIR}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_unique_SO_ID_16_09_25_so_id_mapping_with_assembly_info.tsv" \
 --ref_id_mapping "${MAPDIR}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_with_Uniprot_408_ref_unique_SO_ID_with_assembly_info_ref_id_mapping.tsv" \
 --cell_type "huvec" \
 --outdir "$OUTDIR"


python peptide_groups_pd_analysis_with_reference_19_09_25.py \
 --peptides_file "${DATADIR}/20250505_AS_LC4_MAA_20050_01_VLD_iPSC_F_PeptideGroups.txt" \
 --so_id_mapping_file "${MAPDIR}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_unique_SO_ID_16_09_25_so_id_mapping_with_assembly_info.tsv" \
  --ref_id_mapping "${MAPDIR}/Unique_proteins_Masspec_NMD_RI_HUVEC_CM_tama_with_Uniprot_408_ref_unique_SO_ID_with_assembly_info_ref_id_mapping.tsv" \
 --cell_type "cm" \
 --outdir "$OUTDIR"


# Only need the genomic positions for HUVEC right now as for CM there is no Riboseq data
# Prepare per assembly the unique peptides (no redundancy) to calculate genomic positions 
# with the Split-ORF pipeline
python prepare_unique_peptides_for_genomic_positions.py \
 --unique_peptide_information_csv ${OUTDIR}/huvec_validated_SO_protein_original_Ids_with_assembly.csv \
 --favorite_assembly TAMA_HUVEC \
 --cell_type HUVEC



# Get all .sh files into an array
unique_pep_coord_files=(
    "${OUTDIR}/coordinates/TAMA_HUVEC_unique_peptides_of_HUVEC.bed" 
    "${OUTDIR}/coordinates/TAMA_CM_unique_peptides_of_HUVEC.bed" 
    "${OUTDIR}/coordinates/NMD_unique_peptides_of_HUVEC.bed" 
    "${OUTDIR}/coordinates/RI_unique_peptides_of_HUVEC.bed" 
    )

coord_files=(
    "${SO_PATH}/run_12.09.2025-17.51.04_HUVEC_tama_merged/HUVEC_merged_tama_gene_id_ExonCoordsOfTranscriptsForSO.txt_transcript_positions.bed"
    "${SO_PATH}/run_12.09.2025-14.10.14_CM_tama_merged/CM_merged_tama_gene_id_ExonCoordsOfTranscriptsForSO.txt_transcript_positions.bed"
    "${SO_PATH}/run_30.09.2025-11.30.56_NMD_transcripts_correct_TSL_ref/ExonCoordsWIthChr110_transcript_positions.bed"
    "${SO_PATH}/run_30.09.2025-11.30.56_NMD_transcripts_correct_TSL_ref/ExonCoordsWIthChr110_transcript_positions.bed"
    )



# Iterate through them
for i in "${!unique_pep_coord_files[@]}"; do

    unique_pep_file=${unique_pep_coord_files[$i]}
    coord_file=${coord_files[$i]}
    OUTDIR=$(dirname $unique_pep_file)
    outname=$(basename $unique_pep_file .bed)
    python /home/ckalk/tools/SplitORF_pipeline/Genomic_scripts_18_10_24/genomic_DNA_regions_polars.py \
     $unique_pep_file \
     $coord_file\
    $OUTDIR/${outname}_genomic_coordinates.bed

done

# after the conversion
# python helper_scripts/get_ensembl_gene_ids_from_uniprot_mapping.py \
#  --uniprot_ensembl_tsv ${OUTDIR}/background/Uniprot_to_ensembl_all_genes_HUVEC.tsv

# After GO analysis with SO genes (server)
# python /Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/gprofiler_scripts/create_grprofiler_barplot.py \
#  ${OUTDIR}/validated_genes_for_go/gProfiler_hsapiens_01-10-2025_15-26-41__intersections_high_set_expressed_proteins_genes_background_for_plotting.csv \
#  ${OUTDIR}/validated_genes_for_go/gprofiler_nice_barplot.svg \
# ";"