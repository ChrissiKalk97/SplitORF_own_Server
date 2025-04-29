import sys
import os
import os.path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt


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


def plot_summary(summary_mapping_stats_df, outdir, raw=False):
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

    # rename columns for better readability
    summary_mapping_stats_df_t = summary_mapping_stats_df_t.rename(columns={'Ens_Gencode_lncRNA_ncRNA': 'ncRNA',
                                                                            'Homo_sapiens.GRCh38.dna.chromosome.MT': 'MT',
                                                                            'MANE.GRCh38.v0.95.select_ensembl_rna': 'MANE_mRNA',
                                                                            'hg38-tRNAs': 'tRNA',
                                                                            'rRNA_ref_NCBI_Ens': 'rRNA'})
    # Make one pd.Series per source
    mRNA_df = summary_mapping_stats_df_t["MANE_mRNA"]
    rRNA_df = summary_mapping_stats_df_t["rRNA"]
    ncRNA_df = summary_mapping_stats_df_t["ncRNA"]
    mt_df = summary_mapping_stats_df_t["MT"]
    tRNA_df = summary_mapping_stats_df_t["tRNA"]

    x = ['_'.join(x[2:8]) for x in mRNA_df.index.str.split('_')]

    plt.bar(x, ncRNA_df,
            color='#40A3CD', label='ncRNA')
    plt.bar(x, tRNA_df, bottom=ncRNA_df, color='#019c91', label='tRNA')
    plt.bar(x, rRNA_df, bottom=ncRNA_df + tRNA_df,
            color='#b79165', label='rRNA')
    plt.bar(x, mt_df, bottom=ncRNA_df + tRNA_df + rRNA_df,
            color='orange', label='MT')
    plt.bar(x, mRNA_df, bottom=ncRNA_df + tRNA_df +
            rRNA_df + mt_df, color='#c73832', label='mRNA')

    # add percentage labels
    add_labels(x, mRNA_df, ncRNA_df + tRNA_df + rRNA_df + mt_df)
    add_labels(x, mt_df, ncRNA_df + tRNA_df + rRNA_df)
    add_labels(x, rRNA_df, ncRNA_df + tRNA_df)
    add_labels(x, tRNA_df, ncRNA_df)
    add_labels(x, ncRNA_df, None)

    if raw:
        unmapped_filtered_df = summary_mapping_stats_df_t["filtered_unmapped"]
        plt.bar(x, unmapped_filtered_df, bottom=ncRNA_df +
                tRNA_df + rRNA_df + mt_df + mRNA_df,
                color='#FDDA0D', label='filtered_unmapped')
        add_labels(x, unmapped_filtered_df, ncRNA_df +
                   tRNA_df + rRNA_df + mt_df + mRNA_df)

    plt.ylabel('Percentage')
    plt.xticks(rotation=90)
    plt.legend(title='RNA Type', loc='upper left', bbox_to_anchor=(1, 1))
    plt.title('Stacked Mapping Percentages per Sample')

    if raw:
        plt.savefig(os.path.join(
            outdir, "barplot_mapping_stats_per_sample_raw.png"), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(
            outdir, "barplot_mapping_stats_per_sample.png"), bbox_inches='tight')
    plt.close()


def summarize_mapping_statistics(reference_dict, idx_stats, raw=False):
    """
    Obtain mapping percentages of references from samtools idxstats.

    Each mutlifasta file given as input (list with commaseparated fasta files) is a different source.
    Summarize the alignments by source from mapping stats of single transcripts.

    """
    raw_read_counts_dict = {'uf_muellermcnicoll_2025_04_05_huvec_dhypo_3': 29144203,
                            'uf_muellermcnicoll_2025_04_01_huvec_dnor_2': 25130356,
                            'uf_muellermcnicoll_2025_04_02_huvec_dnor_3': 23952402,
                            'uf_muellermcnicoll_2025_04_03_huvec_dnor_4': 29872996,
                            'uf_muellermcnicoll_2025_04_04_huvec_dhypo_2': 26795853,
                            'uf_muellermcnicoll_2025_04_06_huvec_dhypo_4': 29422604
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

        # if raw:
        #     idx_stat_df_by_source[f'{idx_name}_percentage'] = idx_stat_df_by_source['mapped'] / \
        #         raw_read_counts_dict[idx_name]
        #     idx_stat_df_by_source[f'{idx_name}_not_mapped_or_filtered']\
        #         = raw_read_counts_dict[idx_name] - idx_stat_df_by_source['mapped']
        #     idx_stat_df_by_source[f'{idx_name}_not_mapped_or_filtered_percentage']\
        #         = (raw_read_counts_dict[idx_name] - idx_stat_df_by_source['mapped'])/raw_read_counts_dict[idx_name]
        # else:

        ############### PLOT VIOLIN PLOTS INDIVIDUALLY ##########################################################
        fig, ax = plt.subplots()
        ax.pie(idx_stat_df_by_source['mapped'],
               labels=idx_stat_df_by_source.index, autopct='%1.1f%%')
        plt.title(f'Mapping statistics of {idx_name}')
        if raw:
            plt.savefig(
                os.path.join(os.path.dirname(idx_stat), f'{idx_name}_mapping_percentages_raw.png'))
        else:
            plt.savefig(
                os.path.join(os.path.dirname(idx_stat), f'{idx_name}_mapping_percentages.png'))
        plt.close()

    plot_summary(summary_mapping_stats_df, os.path.dirname(idx_stat), raw)

    return summary_mapping_stats_df


def main():
    ref_fastas = sys.argv[1]
    idx_dir = sys.argv[2]

    idx_stats = []

    for filename in os.listdir(idx_dir):
        if filename.endswith(".out"):
            idx_stats.append(os.path.join(idx_dir, filename))

    reference_dict = assign_trans_IDs_to_source(ref_fastas)

    summarized_map_stats_df = summarize_mapping_statistics(
        reference_dict, idx_stats, raw=True)  #

    out_dir = os.path.dirname(idx_stats[0])
    summarized_map_stats_df.to_csv(f'{out_dir}/summarized_map_stats_df.csv')


if __name__ == '__main__':
    main()
