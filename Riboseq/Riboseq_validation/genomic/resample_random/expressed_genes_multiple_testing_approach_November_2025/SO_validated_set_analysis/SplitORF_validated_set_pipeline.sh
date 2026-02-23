#!/bin/bash -l
# this script runs the downstream analysis steps for the validated regions per sample
# 1. get all validated regions union and statistics for NMD and RI
# 2. run DISGENET analysis for all validated URs (per sample and union/intersection)
# NOTE: DISGENET results first need to be acquired and positioned in the respective folders
# 3. check for RBP validation from RBPDB: the database needs to be downloaded and paths 
# indicated for this to run


eval "$(conda shell.bash hook)"
conda activate Riboseq

disgenet=false
rbpdb=true

outdir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation"
if [[ ! -d "${outdir}/NMD" ]]; then
    python so_validation_pipeline.py \
        --so_results "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-13.07.17_NMD_for_paper/UniqueProteinORFPairs.txt" \
        --ribo_coverage_path "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome" \
        --ur_path "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-13.07.17_NMD_for_paper/Unique_DNA_Regions_genomic_final.bed" \
        --result_dir "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis" \
        --region_type "NMD"
fi

if [[ ! -d "${outdir}/RI" ]]; then
    python so_validation_pipeline.py \
        --so_results "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-14.09.13_RI_for_paper/UniqueProteinORFPairs.txt" \
        --ribo_coverage_path "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/RI_genome" \
        --ur_path "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-14.09.13_RI_for_paper/Unique_DNA_Regions_genomic_final.bed" \
        --result_dir "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis" \
        --region_type "RI"
    fi

if [[ "${disgenet}" == true ]]; then
    bash run_DISGENET_analysis.sh > run_DISGENET_analysis.out 2>&1
fi

if [[ "${rbpdb}" == true ]]; then 
    bash ../downstream_analysis_validated_URs/run_rbpdb_analysis.sh > ../downstream_analysis_validated_URs/outreports_of_runs/run_rbpdb_analysis.out 2>&1
fi

