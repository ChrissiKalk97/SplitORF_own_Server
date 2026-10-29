#!/bin/bash

#----- This script intersects the ultraconserved regions from literature with  ----- #
# ----- Split-ORf unique regions, stop codons and start codons            ----- #

eval "$(conda shell.bash hook)"
conda activate Riboseq

#split_orf_unique_regions="/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_CM_cont_subtraction/Unique_DNA_Regions_genomic_final.bed"

uce_regions="/projects/splitorfs/work/split-orf-prediction/ultraconserved_regions/UCE_regions.bed"
uce_regions_sorted="/projects/splitorfs/work/split-orf-prediction/ultraconserved_regions/UCE_regions_sorted.bed"

ribocov_genes_array=(
    "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome/SO_coverage_categorization/combined_CM_interesting_candidate_genes.txt"
    "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome/SO_coverage_categorization/combined_HUVEC_interesting_candidate_genes.txt"
 )

ribocov_candidate_array=(
    "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome/SO_coverage_categorization/CM_DHYPO/DHYPO_CM_interesting_candidates.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome/SO_coverage_categorization/CM_DNOR/DNOR_CM_interesting_candidates.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome/SO_coverage_categorization/CM_ERR/ERR_CM_interesting_candidates.csv"
    "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome/SO_coverage_categorization/HUVEC_DHYPO/DHYPO_HUVEC_interesting_candidates.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome/SO_coverage_categorization/HUVEC_DNOR/DNOR_HUVEC_interesting_candidates.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome/SO_coverage_categorization/HUVEC_endothel_control/endothel_control_HUVEC_interesting_candidates.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome/SO_coverage_categorization/HUVEC_endothel_treatment/endothel_treatment_HUVEC_interesting_candidates.csv"
)

ur_covered_df_array=(
    "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome/SO_coverage_categorization/CM_DHYPO/val_dna_overlapping_ur_df.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome/SO_coverage_categorization/CM_DNOR/val_dna_overlapping_ur_df.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome/SO_coverage_categorization/CM_ERR/val_dna_overlapping_ur_df.csv"
    "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome/SO_coverage_categorization/HUVEC_DHYPO/val_dna_overlapping_ur_df.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome/SO_coverage_categorization/HUVEC_DNOR/val_dna_overlapping_ur_df.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome/SO_coverage_categorization/HUVEC_endothel_control/val_dna_overlapping_ur_df.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome/SO_coverage_categorization/HUVEC_endothel_treatment/val_dna_overlapping_ur_df.csv"
)

ur_array=(
    "/projects/splitorfs/work/split-orf-prediction/Output/run_07.10.2026-10.43.25_CM_assembly_ref_ORFanage_FIRST_and_Ens110/Unique_DNA_Regions_genomic_final.bed"
    "/projects/splitorfs/work/split-orf-prediction/Output/run_07.10.2026-11.40.34_HUVEC_ref_ORFanage_FIRST_and_Ens110/Unique_DNA_Regions_genomic_final.bed"
 )

for  i in "${!ur_array[@]}"; do
    split_orf_unique_regions=${ur_array[$i]}
    ribocov_genes=${ribocov_genes_array[$i]}
    ribocov_candidates=${ribocov_candidate_array[$i]}
    urs_ribocov=${ur_covered_df_array[$i]}

    outdir=$(dirname "$split_orf_unique_regions")
    sort -k1,1 -k2,2n $uce_regions -o $uce_regions_sorted
    mkdir -p "$outdir"/ultraconservation
    bedtools intersect -a "$split_orf_unique_regions" -b "$uce_regions_sorted" -wao > "$outdir"/ultraconservation/Unique_DNA_Regions_genomic_final_UCE_regions_intersection.bed
    python analyze_ultrconservation_ur.py \
     "$outdir"/ultraconservation/Unique_DNA_Regions_genomic_final_UCE_regions_intersection.bed \
     "$ribocov_genes" \
     "$ribocov_candidates" \
     "$urs_ribocov"

done

