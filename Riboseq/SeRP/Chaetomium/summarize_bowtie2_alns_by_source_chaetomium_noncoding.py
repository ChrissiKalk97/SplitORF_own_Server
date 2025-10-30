import sys
import os
import os.path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def assign_trans_IDs_to_source(ref_fastas):
    """
    Assign the IDs in reference fasta files to the source

    Return a dict with the sources as keys and the transcript IDs as values

    """
    ref_fastas = ref_fastas.split(",")
    reference_dict = {}
    for fasta in ref_fastas:
        ref_seqs = SeqIO.parse(fasta, "fasta")
        source = os.path.basename(fasta)
        source = source.rsplit('.', 1)[0]
        for seq in ref_seqs:
            reference_dict[seq.id] = source
    return reference_dict


def plot_summary(summary_mapping_stats_df, outdir, raw=False, figsize=(12, 6)):
    """
        plot summary statistics of the mapping to Ignolia like reference

        barplot colored by reference type with the sample on the x-axis
    """

    def add_labels(x, y, bottom):
        for i in range(len(x)):
            if bottom is not None:
                plt.text(i, y[i] // 2 + bottom[i], round(y[i], 2), ha='center')
                # Placing text at half the bar height
            else:
                plt.text(i, y[i] // 2, round(y[i], 2), ha='center')

    summary_mapping_stats_df_t = summary_mapping_stats_df.transpose()
    summary_mapping_stats_df_t = summary_mapping_stats_df_t * 100
    summary_mapping_stats_df_t = summary_mapping_stats_df_t[summary_mapping_stats_df_t.index.str.contains(
        'percentage')]
    summary_mapping_stats_df_t = summary_mapping_stats_df_t.sort_index()

    # rename columns for better readability
    summary_mapping_stats_df_t = summary_mapping_stats_df_t.rename(columns={'Chaetomium_thermophilum_protein_coding': 'protein_coding',
                                                                            'Chaetomium_thermophilum_noncoding': 'noncoding'})
    # Make one pd.Series per source
    mRNA_df = summary_mapping_stats_df_t["protein_coding"]
    noncoding_df = summary_mapping_stats_df_t["noncoding"]

    sample_names = ['_'.join(x[2:8]) for x in mRNA_df.index.str.split('_')]

    x = np.arange(len(sample_names))  # Numeric x-axis positions

    plt.figure(figsize=figsize)
    plt.bar(sample_names, noncoding_df,
            color='#40A3CD', label='non-coding')

    plt.bar(x, mRNA_df, bottom=noncoding_df, color='#c73832', label='coding')

    # add percentage labels
    add_labels(x, mRNA_df,  noncoding_df)
    add_labels(x, noncoding_df, None)

    if raw:
        unmapped_filtered_df = summary_mapping_stats_df_t["filtered_unmapped"]
        plt.bar(x, unmapped_filtered_df, bottom=noncoding_df + mRNA_df,
                color='#FDDA0D', label='filtered_unmapped')
        add_labels(x, unmapped_filtered_df,
                   noncoding_df + mRNA_df)

    plt.ylabel('Percentage')
    plt.xticks(x, sample_names, rotation=90)
    plt.legend(title='RNA Type', loc='upper left', bbox_to_anchor=(1, 1))
    plt.title('Stacked Mapping Percentages per Sample')

    if raw:
        plt.savefig(os.path.join(
            outdir, "barplot_mapping_stats_per_sample_raw.png"), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(
            outdir, "barplot_mapping_stats_per_sample.png"), bbox_inches='tight')
    plt.close()


def summarize_mapping_statistics(reference_dict, idx_stats, out_dir, raw=False):
    """
    Obtain mapping percentages of references from samtools idxstats.

    Each mutlifasta file given as input (list with commaseparated fasta files) is a different source.
    Summarize the alignments by source from mapping stats of single transcripts.

    """

    raw_read_counts_dict = {'uf_muellermcnicoll_2025_05_01_MD_In_W1': 41187026,
                            'uf_muellermcnicoll_2025_05_02_MD_In_W2': 37812532,
                            'uf_muellermcnicoll_2025_05_03_MD_In_W3': 35027207,
                            'uf_muellermcnicoll_2025_05_04_MD_In_W4': 37592792,

                            'uf_muellermcnicoll_2025_05_05_MD_In_S1': 31865548,
                            'uf_muellermcnicoll_2025_05_06_MD_In_S2': 36430517,
                            'uf_muellermcnicoll_2025_05_07_MD_In_S3': 28792928,
                            'uf_muellermcnicoll_2025_05_08_MD_In_S4': 33322829,

                            'uf_muellermcnicoll_2025_05_09_MD_E_W1': 34491949,
                            'uf_muellermcnicoll_2025_05_10_MD_E_W2': 35741860,
                            'uf_muellermcnicoll_2025_05_11_MD_E_W3': 40716640,
                            'uf_muellermcnicoll_2025_05_12_MD_E_W4': 33632875,

                            'uf_muellermcnicoll_2025_05_13_MD_E_S1': 22922968,
                            'uf_muellermcnicoll_2025_05_14_MD_E_S2': 22032425,
                            'uf_muellermcnicoll_2025_05_15_MD_E_S3': 37880596,
                            'uf_muellermcnicoll_2025_05_16_MD_E_S4': 30767602,

                            'uf_muellermcnicoll_2025_05_17_MD_In_W1_ND': 30744814
                            }

    summary_mapping_stats_df = None
    for idx_stat in idx_stats:
        idx_stat_df = pd.read_csv(idx_stat,
                                  sep='\t',
                                  index_col=False,
                                  header=None,
                                  names=['tID', 'length', 'mapped', 'unmapped'])
        idx_name = os.path.basename(idx_stat).split('.')[0]
        # skip the last row that is just a star and two dots
        idx_stat_df = idx_stat_df.iloc[:-1, :]
        # map the tids to their source
        idx_stat_df['source'] = idx_stat_df['tID'].map(reference_dict)
        if raw:

            nr_unmapped = raw_read_counts_dict[idx_name] - \
                sum(idx_stat_df['mapped'])
            add_unmapped_row = {'tID': 'unmapped', 'length': 0,
                                'mapped': nr_unmapped, 'unmapped': 0, 'source': 'filtered_unmapped'}
            idx_stat_df = pd.concat(
                [idx_stat_df, pd.DataFrame([add_unmapped_row])], ignore_index=True)

        # get the number of mapping reads per tID
        idx_stat_df_by_source = idx_stat_df[[
            'source', 'mapped']].groupby('source').agg('sum')

        # calculate the mapping percentage per source of all mapping reads
        # if raw equals to true than all the unmapped and filtered reads are counted as well
        idx_stat_df_by_source[f'{idx_name}_percentage'] = idx_stat_df_by_source['mapped'] / \
            sum(idx_stat_df_by_source['mapped'])

        if summary_mapping_stats_df is not None:
            summary_mapping_stats_df[f'{idx_name}_percentage'] = idx_stat_df_by_source[f'{idx_name}_percentage']
            summary_mapping_stats_df[f'{idx_name}_mapped'] = idx_stat_df_by_source['mapped']
        else:
            summary_mapping_stats_df = idx_stat_df_by_source.rename(
                columns={'mapped': f'{idx_name}_mapped'})

        ############### PLOT VIOLIN PLOTS INDIVIDUALLY ##########################################################
        fig, ax = plt.subplots()
        ax.pie(idx_stat_df_by_source['mapped'],
               labels=idx_stat_df_by_source.index, autopct='%1.1f%%')
        plt.title(f'Mapping statistics of {idx_name}')
        if raw:
            plt.savefig(
                os.path.join(out_dir, f'{idx_name}_mapping_percentages_raw.png'))
        else:
            plt.savefig(
                os.path.join(out_dir, f'{idx_name}_mapping_percentages.png'))
        plt.close()

    plot_summary(summary_mapping_stats_df, out_dir, raw)

    return summary_mapping_stats_df


def main():
    ref_fastas = sys.argv[1]
    idx_dir = sys.argv[2]

    idx_stats = []

    for filename in os.listdir(idx_dir):
        if filename.endswith("idxstats.out"):
            idx_stats.append(os.path.join(idx_dir, filename))

    out_dir = os.path.dirname(idx_stats[0])

    reference_dict = assign_trans_IDs_to_source(ref_fastas)

    summarized_map_stats_df = summarize_mapping_statistics(
        reference_dict, idx_stats, out_dir=out_dir)
    summarized_map_stats_df = summarize_mapping_statistics(
        reference_dict, idx_stats, out_dir=out_dir, raw=True)

    summarized_map_stats_df.to_csv(f'{out_dir}/summarized_map_stats_df.csv')

    sys.stdout.flush()
    sys.stderr.flush()


if __name__ == '__main__':
    main()
