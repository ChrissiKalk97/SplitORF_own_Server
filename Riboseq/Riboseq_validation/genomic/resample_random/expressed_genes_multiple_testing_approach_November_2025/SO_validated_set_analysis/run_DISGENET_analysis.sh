#!/bin/bash -l

eval "$(conda shell.bash hook)"
conda activate Riboseq

################################ NMD ###########################################
# glioblastoma against glioblastoma
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/glioblastoma/search_result_C0017636-C1621958-C1514422-C0278878-C0349543+5.xlsx' \
    --samples 'SRR10' \
    --region_type 'NMD' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf' \
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results


# fibrocystic disease
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Breast_Fibrocystic_Disease/search_result_C0016034-C1332629-C0016033-C1971816-C1527396+5.xlsx' \
    --samples 'SRR859076' \
    --region_type 'NMD' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf' \
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results

# fibrocystic disease samples, against malignant neoplasm of breast
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Malignant_neoplasm_of_breast/search_result_C0006142-C1261325-C5441655-C0235653-C0242787+5.xlsx' \
    --samples 'SRR859076' \
    --region_type 'NMD' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf' \
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results



# hmec and T47D only for RI, not found in NMD

# ZR75-1 samples, against malignant neoplasm of breast
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Malignant_neoplasm_of_breast/search_result_C0006142-C1261325-C5441655-C0235653-C0242787+5.xlsx' \
    --samples 'SRR8590790' \
    --region_type 'NMD' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf' \
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results


################################ RI ###########################################
# glioblastoma against glioblastoma
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/glioblastoma/search_result_C0017636-C1621958-C1514422-C0278878-C0349543+5.xlsx' \
    --samples 'SRR10' \
    --region_type 'RI' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf' \
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results


# fibrocystic disease
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Breast_Fibrocystic_Disease/search_result_C0016034-C1332629-C0016033-C1971816-C1527396+5.xlsx' \
    --samples 'SRR859076' \
    --region_type 'RI' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf'\
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results

# fibrocystic disease samples, against malignant neoplasm of breast
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Malignant_neoplasm_of_breast/search_result_C0006142-C1261325-C5441655-C0235653-C0242787+5.xlsx' \
    --samples 'SRR859076' \
    --region_type 'RI' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf'\
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results

# hmec samples, against malignant neoplasm of breast
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Malignant_neoplasm_of_breast/search_result_C0006142-C1261325-C5441655-C0235653-C0242787+5.xlsx' \
    --samples 'SRR859075' \
    --region_type 'RI' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf' \
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results

# T47D (samples kept) samples, against malignant neoplasm of breast
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Malignant_neoplasm_of_breast/search_result_C0006142-C1261325-C5441655-C0235653-C0242787+5.xlsx' \
    --samples 'SRR859077' \
    --region_type 'RI' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf' \
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results

# ZR75-1 samples, against malignant neoplasm of breast
python disgenet_analysis.py \
    --validated_so_csv '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI/validated_so_df.csv' \
    --disgenet_excel '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Malignant_neoplasm_of_breast/search_result_C0006142-C1261325-C5441655-C0235653-C0242787+5.xlsx' \
    --samples 'SRR8590790' \
    --region_type 'RI' \
    --gtf '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf'\
    --result_dir /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/Ribo_seq_intersection_results