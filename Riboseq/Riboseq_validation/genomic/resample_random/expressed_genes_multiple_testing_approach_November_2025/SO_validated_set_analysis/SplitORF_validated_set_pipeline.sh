#!/bin/bash -l

eval "$(conda shell.bash hook)"
conda activate Riboseq

python so_validation_pipeline.py \
    --so_results "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-13.07.17_NMD_for_paper/UniqueProteinORFPairs.txt" \
    --ribo_coverage_path "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome" \
    --ur_path "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-13.07.17_NMD_for_paper/Unique_DNA_Regions_genomic_final.bed" \
    --result_dir "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis" \
    --region_type "NMD"


python so_validation_pipeline.py \
    --so_results "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-14.09.13_RI_for_paper/UniqueProteinORFPairs.txt" \
    --ribo_coverage_path "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome" \
    --ur_path "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-14.09.13_RI_for_paper/Unique_DNA_Regions_genomic_final.bed" \
    --result_dir "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis" \
    --region_type "RI"