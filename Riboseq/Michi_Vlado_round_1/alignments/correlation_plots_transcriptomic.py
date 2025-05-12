# idea: take two replicates (maybe need to be supplied as input the exact files)
# get the idx stats counts per transcript
# log-normalize to the total counts (pseudo count?)
# scatter plot and calculate correlation

import sys
import os
import os.path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn
import math
import scipy as sp


def read_idx_stats(idx_stat):
    idx_stat_df = pd.read_csv(idx_stat,
                              sep='\t',
                              index_col=False,
                              header=None,
                              names=['tID', 'length', 'mapped', 'unmapped'])
    idx_name = os.path.basename(idx_stat).split('.')[0]
    # skip the last row that is just a star and two dots
    idx_stat_df = idx_stat_df.iloc[:-1, :]
    # percentage and add a pseudocount
    idx_stat_df[f'{idx_name}_mapping_perc'] = idx_stat_df['mapped'] / \
        sum(idx_stat_df['mapped'])
    idx_stat_df[f'{idx_name}_log_mapping_perc'] = idx_stat_df[f'{idx_name}_mapping_perc'].apply(
        lambda x: math.log(x+0.000001))
    idx_stat_df[f'{idx_name}_log_mapping_reads'] = idx_stat_df['mapped'].apply(
        lambda x: math.log(x+0.000001))
    idx_stat_df = idx_stat_df[[
        'tID', f'{idx_name}_mapping_perc', f'{idx_name}_log_mapping_reads', f'{idx_name}_log_mapping_perc']]
    return idx_stat_df, idx_name


def main():
    idx_rep1 = sys.argv[1]
    idx_rep2 = sys.argv[2]
    out_dir = sys.argv[3]

    idx_stat_df_1, idx_name1 = read_idx_stats(idx_rep1)
    idx_stat_df_2, idx_name2 = read_idx_stats(idx_rep2)

    idx_merged_df = pd.merge(
        idx_stat_df_1, idx_stat_df_2, on='tID', how='outer')

    # calculate pearson correlation
    r, _ = sp.stats.pearsonr(
        x=idx_merged_df[f'{idx_name1}_log_mapping_perc'], y=idx_merged_df[f'{idx_name2}_log_mapping_perc'])
    scat_plot = sbn.scatterplot(data=idx_merged_df,
                                x=f'{idx_name1}_log_mapping_perc', y=f'{idx_name2}_log_mapping_perc')

    scat_plot.set(xlabel=f"{idx_name1} log percentage", ylabel=f"{idx_name2} log percentage",
                  title='Scatterplot of log scaled mapping percentage of genes')
    second_rep = idx_name2.split('_')[-1]
    ax = plt.gca()  # Get a matplotlib's axes instance
    plt.text(.05, .8, "Pearson's r ={:.2f}".format(r), transform=ax.transAxes)
    plt.savefig(os.path.join(
        out_dir, f'scatterplot_{idx_name1}_{second_rep}.png'), bbox_inches='tight', dpi=300)
    plt.close()

    # how to get both in the same df?
    sys.stdout.flush()
    sys.stderr.flush()


if __name__ == '__main__':
    main()
