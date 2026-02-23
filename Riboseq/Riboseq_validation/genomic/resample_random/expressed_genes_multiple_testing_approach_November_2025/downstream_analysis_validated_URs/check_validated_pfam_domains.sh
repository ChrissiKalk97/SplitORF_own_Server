#!/bin/bash -l

eval "$(conda shell.bash hook)"
conda activate Riboseq

################################ NMD ###########################################
python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/expressed_genes_multiple_testing_approach_November_2025/downstream_analysis_validated_URs/check_validated_pfam_domains.py \
    --annotated_pfam_file '/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-13.07.17_NMD_for_paper/UniqueProteinORFPairs_annotated.txt' \
    --validation_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/validated_so_df.csv' \
    --region_type 'NMD' \
    --rbp_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/proteins.php'


################################ RI ###########################################
python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/expressed_genes_multiple_testing_approach_November_2025/downstream_analysis_validated_URs/check_validated_pfam_domains.py \
    --annotated_pfam_file '/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-14.09.13_RI_for_paper/UniqueProteinORFPairs_annotated.txt' \
    --validation_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI/validated_so_df.csv' \
    --region_type 'RI' \
    --rbp_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/proteins.php'
