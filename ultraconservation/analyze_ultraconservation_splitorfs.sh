#!/bin/bash

#----- This script intersects the ultraconserved regions from literature with  ----- #
# ----- Split-ORf unique regions, stop codons and start codons            ----- #

eval "$(conda shell.bash hook)"
conda activate Riboseq

#split_orf_unique_regions="/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed"

uce_regions="/projects/splitorfs/work/split-orf-prediction/ultraconserved_regions/UCE_regions.bed"
uce_regions_sorted="/projects/splitorfs/work/split-orf-prediction/ultraconserved_regions/UCE_regions_sorted.bed"

ribocov_genes_array=("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/combined_NMD_interesting_candidate_genes.txt" "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/RI_genome/SO_coverage_categorization/combined_RI_interesting_candidate_genes.txt")
ur_array=("/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed" "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/Unique_DNA_Regions_genomic_final.bed")
for  i in "${!ur_array[@]}"; do
    split_orf_unique_regions=${ur_array[$i]}
    ribocov_genes=${ribocov_genes_array[$i]}
    outdir=$(dirname "$split_orf_unique_regions")
    sort -k1,1 -k2,2n $uce_regions -o $uce_regions_sorted
    mkdir -p "$outdir"/ultraconservation
    bedtools intersect -a "$split_orf_unique_regions" -b "$uce_regions_sorted" -wao > "$outdir"/ultraconservation/Unique_DNA_Regions_genomic_final_UCE_regions_intersection.bed
    python analyze_ultrconservation_ur.py \
     "$outdir"/ultraconservation/Unique_DNA_Regions_genomic_final_UCE_regions_intersection.bed \
     "$ribocov_genes"
done

