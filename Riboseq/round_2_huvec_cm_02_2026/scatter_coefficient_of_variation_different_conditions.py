# calculate coefficient of variation for each transcript
# average across samples, report median CV and plot a histogram
import os
import argparse

import pyBigWig
import numpy as np
import pandas as pd

import seaborn as sbn
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('--path_to_cv_csvs_1',
                        help='Path to CSV files with the CV values per transcript of round 1')
    parser.add_argument('--path_to_cv_csvs_2',
                        help='Path to CSV files with the CV values per transcript of round 2')
    parser.add_argument('--out_path', help='Path where to output the plots')
    parser.add_argument('--condition1', help='condition to be plotted')
    parser.add_argument('--condition2', help='condition to be plotted')
    parser.add_argument('--color', help='color of plot')

    return parser.parse_args()


# path_to_cv_csvs_1 = '/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files/trans_cv_dfs'
# path_to_cv_csvs_2 = '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files/trans_cv_dfs'
# transcript_fai = '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai'
# mane_transcripts_cds_bed = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed'
# condition = 'DNOR'
# out_path = '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files/metagene_plots'
# color = '#1eb0e6'


def main(path_to_cv_csvs_1, path_to_cv_csvs_2,
         out_path, condition1, condition2, color='blue'):
    def plot_scatter(cv_df_concat, out_path, condition1, condition2, color):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=600)
        # plt.figure(figsize=(10, 5))
        sbn.scatterplot(cv_df_concat, x='cv round1',
                        y='cv round2', size=0.5, legend=False)

        ax.plot(
            *np.array([ax.get_xlim(), np.array(ax.get_xlim())]), color="red")

        ax.set_title(
            f'coefficient of variation in transcripts in {condition2} vs {condition1}')

        plt.xlim(0, 10)
        plt.ylim(0, 10)

        plt.savefig(os.path.join(out_path, f'{condition2}_{condition1}_cv_scatter.svg'),
                    format='svg', dpi=300, bbox_inches='tight')
        plt.close()

    for csv in os.listdir(path_to_cv_csvs_1):
        if csv.endswith('.csv') and condition1 in csv.upper():
            cv_df_1 = pd.read_csv(os.path.join(path_to_cv_csvs_1, csv))
            cv_df_1 = cv_df_1[cv_df_1['cv'] > 0]

    for csv in os.listdir(path_to_cv_csvs_2):
        if csv.endswith('.csv') and condition2 in csv.upper():
            cv_df_2 = pd.read_csv(os.path.join(path_to_cv_csvs_2, csv))
            cv_df_2 = cv_df_2[cv_df_2['cv'] > 0]

    tids_keep = [tid for tid in cv_df_1['tid']
                 if tid in cv_df_2['tid'].tolist()]
    print(f'Nr of transcripts with CPM > 5 in both datasets {condition2} and {condition1}:', len(
        tids_keep))

    cv_df_1 = cv_df_1[cv_df_1['tid'].isin(tids_keep)].copy()
    cv_df_2 = cv_df_2[cv_df_2['tid'].isin(tids_keep)].copy()

    cv_df_1 = cv_df_1.set_index('tid')
    cv_df_2 = cv_df_2.set_index('tid')

    cv_df_1 = cv_df_1.rename(columns={'cv': 'cv round1'})
    cv_df_2 = cv_df_2.rename(columns={'cv': 'cv round2'})

    cv_df_concat = pd.concat(
        [cv_df_1['cv round1'], cv_df_2['cv round2']], axis=1)

    plot_scatter(cv_df_concat, out_path, condition1, condition2, color)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_cv_csvs_1 = args.path_to_cv_csvs_1
    path_to_cv_csvs_2 = args.path_to_cv_csvs_2
    out_path = args.out_path
    condition1 = args.condition1
    condition2 = args.condition2
    color = args.color

    main(path_to_cv_csvs_1, path_to_cv_csvs_2,
         out_path, condition1, condition2, color)
