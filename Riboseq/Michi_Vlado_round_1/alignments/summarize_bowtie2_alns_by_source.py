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

    sample_names = ['_'.join(x[2:9]) if 'RR' in x else '_'.join(
        x[2:8]) for x in mRNA_df.index.str.split('_')]
    plt.figure(figsize=figsize)
    plt.bar(sample_names, ncRNA_df,
            color='#40A3CD', label='ncRNA')
    plt.bar(sample_names, tRNA_df, bottom=ncRNA_df,
            color='#019c91', label='tRNA')
    plt.bar(sample_names, rRNA_df, bottom=ncRNA_df + tRNA_df,
            color='#b79165', label='rRNA')
    plt.bar(sample_names, mt_df, bottom=ncRNA_df + tRNA_df + rRNA_df,
            color='orange', label='MT')
    plt.bar(sample_names, mRNA_df, bottom=ncRNA_df + tRNA_df +
            rRNA_df + mt_df, color='#c73832', label='mRNA')

    # add percentage labels
    add_labels(sample_names, mRNA_df, ncRNA_df + tRNA_df + rRNA_df + mt_df)
    add_labels(sample_names, mt_df, ncRNA_df + tRNA_df + rRNA_df)
    add_labels(sample_names, rRNA_df, ncRNA_df + tRNA_df)
    add_labels(sample_names, tRNA_df, ncRNA_df)
    add_labels(sample_names, ncRNA_df, None)

    if raw:
        unmapped_filtered_df = summary_mapping_stats_df_t["filtered_unmapped"]
        plt.bar(sample_names, unmapped_filtered_df, bottom=ncRNA_df +
                tRNA_df + rRNA_df + mt_df + mRNA_df,
                color='#FDDA0D', label='filtered_unmapped')
        add_labels(sample_names, unmapped_filtered_df, ncRNA_df +
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


def summarize_mapping_statistics(reference_dict, idx_stats, out_dir, raw=False):
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
                            'uf_muellermcnicoll_2025_04_06_huvec_dhypo_4': 29422604,
                            'uf_muellermcnicoll_2025_04_07_In_Puro_1': 24221259,
                            'uf_muellermcnicoll_2025_04_08_In_Puro_3': 29235694,
                            'uf_muellermcnicoll_2025_04_09_In_Puro_4': 26411807,
                            'uf_muellermcnicoll_2025_04_10_In_CHX_1': 27106673,
                            'uf_muellermcnicoll_2025_04_11_In_CHX_2': 27709830,
                            'uf_muellermcnicoll_2025_04_12_In_CHX_4': 26852660,
                            'uf_muellermcnicoll_2025_04_13_IP_Puro_1': 23927046,
                            'uf_muellermcnicoll_2025_04_14_IP_Puro_3': 25064225,
                            'uf_muellermcnicoll_2025_04_15_IP_Puro_4': 24337638,
                            'uf_muellermcnicoll_2025_04_16_IP_CHX_1': 24226505,
                            'uf_muellermcnicoll_2025_04_17_IP_CHX_2': 20622935,
                            'uf_muellermcnicoll_2025_04_18_IP_CHX_4': 28061501,
                            'uf_muellermcnicoll_2025_05_18_RR_A2_CHX_E1': 26281252,
                            'uf_muellermcnicoll_2025_05_19_RR_A2_CHX_E2': 21850590,
                            'uf_muellermcnicoll_2025_05_20_RR_A2_CHX_E3': 28072070,
                            'uf_muellermcnicoll_2025_05_21_RR_A2_CHX_E4': 25134187,
                            'uf_muellermcnicoll_2025_05_22_RR_A2_Puro_E1': 26937848,
                            'uf_muellermcnicoll_2025_05_23_RR_A2_Puro_E2': 21151387,
                            'uf_muellermcnicoll_2025_05_24_RR_A2_Puro_E3': 26771687,
                            'uf_muellermcnicoll_2025_05_25_RR_B1_CHX_E1': 29423782,
                            'uf_muellermcnicoll_2025_05_26_RR_B1_CHX_E2': 28519436,
                            'uf_muellermcnicoll_2025_05_27_RR_B1_CHX_E3': 20533262,
                            'uf_muellermcnicoll_2025_05_28_RR_B1_CHX_E4': 22720343,
                            'uf_muellermcnicoll_2025_05_29_RR_B1_Puro_E1': 22952596,
                            'uf_muellermcnicoll_2025_05_30_RR_B1_Puro_E2': 16497721,
                            'uf_muellermcnicoll_2025_05_31_RR_B1_Puro_E3': 22143053,
                            'uf_muellermcnicoll_2025_05_32_RR_In_CHX_1': 27157334,
                            'uf_muellermcnicoll_2025_05_33_RR_In_CHX_2': 31787851,
                            'uf_muellermcnicoll_2025_05_34_RR_In_CHX_3': 27470164,
                            'uf_muellermcnicoll_2025_05_35_RR_In_CHX_4': 18438108,
                            'uf_muellermcnicoll_2025_05_36_RR_In_Puro_1': 27408841,
                            'uf_muellermcnicoll_2025_05_37_RR_In_Puro_2': 25715829,
                            'uf_muellermcnicoll_2025_05_38_RR_In_Puro_3': 26110767,
                            'uf_muellermcnicoll_2025_05_39_RR_M_CHX_E1': 16836171,
                            'uf_muellermcnicoll_2025_05_40_RR_M_CHX_E2': 21822278,
                            'uf_muellermcnicoll_2025_05_41_RR_M_CHX_E3': 12929426,
                            'uf_muellermcnicoll_2025_05_42_RR_M_CHX_E4': 30557025
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
