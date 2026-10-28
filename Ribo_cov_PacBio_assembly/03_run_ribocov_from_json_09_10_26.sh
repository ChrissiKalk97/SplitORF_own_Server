#!/bin/bash


# =============================================================================
# Script Name: 
# Description: Run the ribocov module on HUVEC and CM in the follwoing sets:
#              - Set 1: CM - DNOR
#              - Set 2: CM - DHYPO
#              - Step 3: 
#              - Step 4: 
#              - Step 5: 
# Usage:       bash 
# Author:      Christina Wrobel
# Date:        2026-10-09
# =============================================================================



eval "$(conda shell.bash hook)"
conda activate split-orf-pipeline-oct-26

################################################################################
# RUN RIBO_COV                                                                 #
################################################################################
# ribo-cov /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/ribo_cov_CM_Ribo_round2_SO_pred_ens_and_ORFanage_FIRST_ref_10_07_26_DNOR.json
# ribo-cov /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/ribo_cov_CM_Ribo_round2_SO_pred_ens_and_ORFanage_FIRST_ref_10_07_26_DHYPO.json
# ribo-cov /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/ribo_cov_HUVEC_Ribo_round2_SO_pred_ens_and_ORFanage_FIRST_ref_10_07_26_DHYPO.json
# ribo-cov /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/ribo_cov_HUVEC_Ribo_round2_SO_pred_ens_and_ORFanage_FIRST_ref_10_07_26_DNOR.json


# preprocessing of these files see here: /home/ckalk/scripts/SplitORFs/Riboseq/New_datasets_riboseq_so_paper/run_Riboseq_analysis_cancer_Riboseq_org.sh
# ribo-cov /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/2026_10_23_ribo_cov_Mauro_data_HUVEC_control.json
# ribo-cov /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/2026_10_23_ribo_cov_Mauro_data_HUVEC_treatment.json
# ribo-cov /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/2026_10_23_ribo_cov_CM_IPSC_public_control.json


################################################################################
# COMBINE THE INTERESTING CANDIDATE GENES WITHIN CELL TYPES                    #
################################################################################
cm_ribocov_dir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome"
huvec_ribocov_dir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome"
find "${cm_ribocov_dir}/SO_coverage_categorization" -type f -name "*genes_interesting_candidates.txt" -print0 | xargs -0 cat | sort | uniq > "${cm_ribocov_dir}/SO_coverage_categorization/combined_CM_interesting_candidate_genes.txt"
find "${huvec_ribocov_dir}/SO_coverage_categorization" -type f -name "*genes_interesting_candidates.txt" -print0 | xargs -0 cat | sort | uniq > "${huvec_ribocov_dir}/SO_coverage_categorization/combined_HUVEC_interesting_candidate_genes.txt"