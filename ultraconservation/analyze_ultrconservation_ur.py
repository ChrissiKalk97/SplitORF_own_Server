# create BED file from the Excel of publicaiton

# conda: SplitORF

import sys
import os.path
import pandas as pd


uce_intersection = sys.argv[1]
ribocov_genes = sys.argv[2]
# uce_intersection = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/ultraconservation/Unique_DNA_Regions_genomic_final_UCE_regions_intersection.bed'
# # uce_intersection_ri = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/ultraconservation/Unique_DNA_Regions_genomic_final_UCE_regions_intersection.bed'
# ribocov_genes = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/combined_NMD_interesting_candidate_genes.txt"
outdir = os.path.dirname(uce_intersection)

uce_intersection_df = pd.read_csv(uce_intersection, sep='\t', header=None)
uce_intersection_df = uce_intersection_df[uce_intersection_df.iloc[:, 10] > 0].reset_index(
    drop=True)

print("Number of intersections with UCEs:", len(uce_intersection_df.index))
# 23 intersections

uce_intersection_df['gene_id'] = uce_intersection_df.iloc[:, 3].apply(
    lambda x: x.split('|')[0])
uce_intersection_df['transcript_id'] = uce_intersection_df.iloc[:, 3].apply(
    lambda x: x.split('|')[1].split(':')[0])

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

print(uce_intersection_df_ribocov)
