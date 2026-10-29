# create BED file from the Excel of publicaiton

# conda: SplitORF

import sys
import os.path
import pandas as pd
import ast

uce_intersection = sys.argv[1]
ribocov_genes = sys.argv[2]
ribocov_candidates_csvs = sys.argv[3]
ribocov_regions_covered_csvs = sys.argv[4]


# uce_intersection = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/ultraconservation/Unique_DNA_Regions_genomic_final_UCE_regions_intersection.bed'
# # uce_intersection_ri = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/ultraconservation/Unique_DNA_Regions_genomic_final_UCE_regions_intersection.bed'
# ribocov_genes = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/combined_NMD_interesting_candidate_genes.txt"

# ribocov_candidates_csvs = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_cancer/cancer_NMD_interesting_candidates.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_HCT_control/HCT_control_NMD_interesting_candidates.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_NMD_inhibition/NMD_inhibition_NMD_interesting_candidates.csv"
# ribocov_regions_covered_csvs = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_cancer/val_dna_overlapping_ur_df.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_HCT_control/val_dna_overlapping_ur_df.csv,/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_NMD_inhibition/val_dna_overlapping_ur_df.csv"

ribocov_candidates_csv_list = ribocov_candidates_csvs.split(",")
ribocov_regions_covered_csv_list = ribocov_regions_covered_csvs.split(",")

region_data_types = [os.path.basename(os.path.dirname(
    csv)) for csv in ribocov_regions_covered_csv_list]

outdir = os.path.dirname(uce_intersection)

uce_intersection_df = pd.read_csv(uce_intersection, sep='\t', header=None, names=[
                                  'chr_UR', 'start_UR', 'stop_UR', 'UCE_overlapping_ORF_ID_for_merge', 'score', 'strand_UR', 'chr_UCE', 'start_UCE', 'stop_UCE', 'UCE_ID', 'bp_overlap'])


uce_intersection_df = uce_intersection_df[uce_intersection_df['bp_overlap'] > 0].reset_index(
    drop=True)

print("Number of intersections with UCEs:", len(uce_intersection_df.index))
# 23 intersections

uce_intersection_df['gene_id'] = uce_intersection_df['UCE_overlapping_ORF_ID_for_merge'].apply(
    lambda x: x.split('|')[0])
uce_intersection_df['transcript_id'] = uce_intersection_df['UCE_overlapping_ORF_ID_for_merge'].apply(
    lambda x: x.split('|')[1].split(':')[0])
uce_intersection_df['orf_id'] = uce_intersection_df['UCE_overlapping_ORF_ID_for_merge'].apply(
    lambda x: x.split(":")[1])

print("Number of UCE intersecting Split-ORF genes:",
      len(uce_intersection_df['gene_id'].unique()))

uce_intersection_df.to_csv(os.path.join(
    outdir, 'uce_intersections.csv'), sep='\t')

ribocov_genes = pd.read_csv(ribocov_genes, header=None)

uce_intersection_df_ribocov = uce_intersection_df[uce_intersection_df['gene_id'].isin(
    ribocov_genes.iloc[:, 0].to_list())].reset_index(drop=True)

print("Number of ribo-cov intersections with UCEs:",
      len(uce_intersection_df_ribocov.index))
print("Number of UCE intersecting ribo-cov Split-ORF genes:",
      len(uce_intersection_df_ribocov['gene_id'].unique()))

uce_intersection_df_ribocov.to_csv(os.path.join(
    outdir, 'uce_intersections_ribocov_genes.csv'), sep='\t')

################################################################################
# Intersection with ribo-cov regions, exact coordinates
################################################################################

uce_intersection_df = uce_intersection_df.rename(
    columns={3: 'UCE_overlapping_ORF_ID_for_merge'})

merged_uce_ribocov_dfs = []

for ribocov_candidate_csv, ribocov_regions_covered_csv, region_data_type in zip(
        ribocov_candidates_csv_list,
        ribocov_regions_covered_csv_list,
        region_data_types):

    ribocov_candidate_covered_df = pd.read_csv(ribocov_candidate_csv)
    ribocov_regions_covered_df = pd.read_csv(
        ribocov_regions_covered_csv, index_col=0)

    ribocov_regions_covered_df['OrfTransID'] = ribocov_regions_covered_df['OrfTransID'].apply(
        ast.literal_eval)

    # relevant for NMD inhibition, subset to ribocov regions of interesting
    # candidates
    so_trans = set(ribocov_candidate_covered_df['OrfTransID'])
    ribocov_regions_covered_candidates_df = ribocov_regions_covered_df[ribocov_regions_covered_df['OrfTransID'].apply(
        lambda x: any(so_id.split('|')[1] in so_trans for so_id in x))].copy()
    ribocov_regions_covered_candidates_df['ID_list'] = ribocov_regions_covered_candidates_df['ID'].apply(
        lambda x: x.split(','))

    uce_urs = set(uce_intersection_df['UCE_overlapping_ORF_ID_for_merge'])
    ribocov_regions_covered_candidates_df['UCE_overlap_ORF'] = ribocov_regions_covered_candidates_df['ID_list'].apply(
        lambda x: [so_id for so_id in x if so_id in uce_urs])
    ribocov_regions_covered_candidates_df['UCE_overlapping_ORF_ID_for_merge'] = ribocov_regions_covered_candidates_df['UCE_overlap_ORF'].apply(
        lambda x: x[0] if len(x) > 0 else None)

    uces_ur_ribocov = ribocov_regions_covered_candidates_df[ribocov_regions_covered_candidates_df['ID_list'].apply(
        lambda x: any(so_id in uce_urs for so_id in x))].copy()

    uce_ribocov_regions_merged = uces_ur_ribocov.merge(
        uce_intersection_df, on='UCE_overlapping_ORF_ID_for_merge', how='inner')

    uce_ribocov_regions_merged['region_type'] = region_data_type

    print(f'Number of distinct URs overlapping UCEs in {region_data_type}:', len(
        uce_ribocov_regions_merged.index))
    print(f'Number of different genes with URs overlapping UCEs in {region_data_type}:', len(
        uce_ribocov_regions_merged['gene_id'].unique()))

    uce_ribocov_regions_merged.to_csv(os.path.join(
        outdir, f'{region_data_type}_uce_ribocov_regions_merged.csv'))

    merged_uce_ribocov_dfs.append(uce_ribocov_regions_merged)


all_merged_df = pd.concat(merged_uce_ribocov_dfs, ignore_index=True)
region_type = region_data_types[0].split('_')[0]

all_merged_df.to_csv(os.path.join(
    outdir, f'{region_type}_uce_ribocov_regions_merged_all.csv'), index=False)

pd.Series(all_merged_df['gene_id'].unique()).to_csv(os.path.join(
    outdir, f'{region_type}_uce_ribocov_regions_genes.txt'), index=False)
