import sys
import os
import os.path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


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


def plot_summary(summary_mapping_stats_df, outdir):
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
    mRNA_df = summary_mapping_stats_df_t["MANE_mRNA"]
    rRNA_df = summary_mapping_stats_df_t["rRNA"]
    ncRNA_df = summary_mapping_stats_df_t["ncRNA"]
    mt_df = summary_mapping_stats_df_t["MT"]
    tRNA_df = summary_mapping_stats_df_t["tRNA"]

    x = ['_'.join(x[2:7]) for x in mRNA_df.index.str.split('_')]

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

    plt.ylabel('Percentage')
    plt.xticks(rotation=45)
    plt.legend(title='RNA Type')
    plt.title('Stacked Mapping Percentages per Sample')

    plt.savefig(os.path.join(outdir, "barplot_mapping_stats_per_sample.png"))
    plt.close()


def summarize_mapping_statistics(reference_dict, idx_stats):
    """
    Obtain mapping percentages of references from samtools idxstats.

    Each mutlifasta file given as input (list with commaseparated fasta files) is a different source.
    Summarize the alignments by source from mapping stats of single transcripts.

    """
    summary_mapping_stats_df = None
    for idx_stat in idx_stats:
        idx_stat_df = pd.read_csv(idx_stat,
                                  sep='\t',
                                  index_col=False,
                                  header=None,
                                  names=['tID', 'length', 'mapped', 'unmapped'])
        idx_name = os.path.basename(idx_stat).split('.')[0]
        idx_stat_df = idx_stat_df.iloc[:-1, :]
        idx_stat_df['source'] = idx_stat_df['tID'].map(reference_dict)
        idx_stat_df_by_source = idx_stat_df[[
            'source', 'mapped']].groupby('source').agg('sum')
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
        plt.savefig(
            os.path.join(os.path.dirname(idx_stat), f'{idx_name}_mapping_percentages.png'))
        plt.close()

    plot_summary(summary_mapping_stats_df, os.path.dirname(idx_stat))

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
        reference_dict, idx_stats)

    out_dir = os.path.dirname(idx_stats[0])
    summarized_map_stats_df.to_csv(f'{out_dir}/summarized_map_stats_df.csv')


if __name__ == '__main__':
    main()
