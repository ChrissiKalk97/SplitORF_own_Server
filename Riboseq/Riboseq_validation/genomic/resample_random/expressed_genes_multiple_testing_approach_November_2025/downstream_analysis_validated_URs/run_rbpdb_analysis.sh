#!/bin/bash -l

eval "$(conda shell.bash hook)"
conda activate Riboseq

################################ NMD ###########################################
python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/expressed_genes_multiple_testing_approach_November_2025/downstream_analysis_validated_URs/rbp_anlaysis_rbpdb.py \
    --rbp_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/proteins.php' \
    --validation_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/validated_so_df.csv' \
    --region_type 'NMD'


################################ RI ###########################################
python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/expressed_genes_multiple_testing_approach_November_2025/downstream_analysis_validated_URs/rbp_anlaysis_rbpdb.py \
    --rbp_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/proteins.php' \
    --validation_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI/validated_so_df.csv' \
    --region_type 'RI'


python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/expressed_genes_multiple_testing_approach_November_2025/downstream_analysis_validated_URs/ri_nmd_union_intersection_rbps.py \
    --ri_validated_rbp_file "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/NMD/rbp_val_df.csv" \
    --nmd_validated_rbp_file "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/RI/rbp_val_df.csv"
