import re
import pandas as pd
import glob
from collections import Counter
import os
import os.path

from plotting import plot_three_category_pie


def validated_so_per_sample_analysis(ribo_coverage_path, all_predicted_so_orfs):
    '''
    for coverage based validated SOs perform background analysis
    '''
    for empirical_Ribo_findings_file in glob.glob(f"{ribo_coverage_path}/*_unique_regions.csv"):
        if 'HCT' in empirical_Ribo_findings_file or \
            'SRR85' in empirical_Ribo_findings_file or \
                'SRR10' in empirical_Ribo_findings_file:
            sample = os.path.basename(
                empirical_Ribo_findings_file).removesuffix('_unique_regions.csv')
            Ribo_results_df = pd.read_csv(
                empirical_Ribo_findings_file, header=0, index_col=0)
            Ribo_results_significant_df = Ribo_results_df[Ribo_results_df["significant"] == 1].copy(
            )
            Ribo_results_significant_df['ORF'] = Ribo_results_significant_df['name'].apply(
                lambda x: x.split(':')[1])
            Ribo_results_significant_df['OrfTransID'] = Ribo_results_significant_df['name'].apply(
                lambda x: x.split(':')[0].split('|')[1])
            all_predicted_so_orfs[sample] = 0
            all_predicted_so_orfs[sample] = all_predicted_so_orfs['OrfID'].isin(
                Ribo_results_significant_df['ORF'])
    return all_predicted_so_orfs


def explode_so_df(predicted_so_orfs):
    predicted_so_orfs = predicted_so_orfs[[
        'OrfTransID', 'OrfIDs', 'OrfStarts', 'geneID']].copy()
    predicted_so_orfs['OrfID'] = predicted_so_orfs.apply(
        lambda x: x['OrfIDs'].split(','), axis=1)
    predicted_so_orfs['OrfStart'] = predicted_so_orfs['OrfStarts']
    all_predicted_so_orfs = predicted_so_orfs.explode(
        ['OrfID', 'OrfStart'], ignore_index=True).copy()
    total_nr_so = len(all_predicted_so_orfs.index)
    total_nr_so_transcripts = len(
        all_predicted_so_orfs['OrfTransID'].unique())
    print('Total number of predicted Split-ORFs', total_nr_so)
    print('Total number of predicted Split-ORF transcripts',
          total_nr_so_transcripts)
    return all_predicted_so_orfs, predicted_so_orfs, total_nr_so


def subset_validated_sos_df(all_predicted_so_orfs, outdir, region_type):
    all_predicted_so_orfs['ValidationCount'] = all_predicted_so_orfs.iloc[:, range(
        6, len(all_predicted_so_orfs.columns))].sum(axis=1, numeric_only=True)
    validated_so_df = all_predicted_so_orfs[all_predicted_so_orfs['ValidationCount'] > 0].copy(
    )

    # write txt file with genes that have validated SOs
    pd.Series(validated_so_df['geneID'].unique()).to_csv(os.path.join(
        outdir, f'SO_validated_genes_{region_type}.txt'), header=False, index=False)
    nr_validated_so = len(validated_so_df.index)
    nr_validated_transcripts = len(validated_so_df['OrfTransID'].unique())
    print('Number of validated SOs:', nr_validated_so)
    print('Number of validated SO transcripts:', nr_validated_transcripts)
    return validated_so_df, nr_validated_so, nr_validated_transcripts


def subset_UR_for_expressed_genes(dna_ur_df, validated_so_df, ribo_coverage_path, outdir, region_type):
    '''
        subset the URs for the genes that are included by TPM filter
    '''
    # get gene and transcript ID in dna_ur_df from name
    dna_ur_df['gene'] = dna_ur_df['OrfTransID'].apply(
        lambda x: x.split('|')[0])
    dna_ur_df['tID'] = dna_ur_df['OrfTransID'].apply(
        lambda x: x.split('|')[1])

    genes_to_keep = filter_so_genes(ribo_coverage_path)
    pd.Series(genes_to_keep).to_csv(os.path.join(
        outdir, f'all_considered_SO_genes_{region_type}.txt'))
    # subset and count ORFs with URs for expressed genes
    dna_ur_df_subset = dna_ur_df[(dna_ur_df['gene'].isin(genes_to_keep)) | (
        dna_ur_df['tID'].isin(validated_so_df['OrfTransID']))]
    nr_orfs_with_UR_expressed = len(dna_ur_df_subset['OrfID'].unique())

    nr_val_SO_trans = len(validated_so_df['OrfTransID'])
    nr_val_SO_trans_filtered = len(
        validated_so_df[validated_so_df['geneID'].isin(genes_to_keep)].index)
    nr_val_lost = nr_val_SO_trans - nr_val_SO_trans_filtered
    # there may not be URs validated of genes not included!
    assert nr_val_lost == 0
    return nr_orfs_with_UR_expressed, dna_ur_df_subset


def filter_so_genes(ribo_coverage_path):
    '''
        filter genes for TPm of at least 20 in at least 2 samples
    '''
    genes_above_20_list = []
    for empirical_Ribo_findings_file in glob.glob(f"{ribo_coverage_path}/**/Unique_DNA_Regions_genomic_*_chrom_sorted.bed", recursive=True):
        if 'HCT' in empirical_Ribo_findings_file or \
            'SRR85' in empirical_Ribo_findings_file or \
                'SRR10' in empirical_Ribo_findings_file:
            unique_regions_tpm_20_filtered = pd.read_csv(
                empirical_Ribo_findings_file, header=None, index_col=None, sep='\t')
            unique_regions_tpm_20_filtered['gene_id'] = unique_regions_tpm_20_filtered.iloc[:, 3].apply(
                lambda x: x.split('|')[0])
            genes_above_20_list = genes_above_20_list + \
                unique_regions_tpm_20_filtered['gene_id'].to_list()
    counts_above_20 = Counter(genes_above_20_list)
    genes_to_keep = []
    for gene, count in counts_above_20.items():
        if count > 0:
            genes_to_keep.append(gene)
    return genes_to_keep


def get_so_position_in_transcript(so_df):
    # sort the ORF starts by position
    so_df['OrfStarts'] = so_df.apply(
        lambda x: sorted([int(start) for start in x['OrfStarts']]), axis=1)
    # map the ORF start to the respective position in the sorted list
    # indicate whether it is the first or a later (first, middle, last)
    so_df['OrfIndex'] = so_df.apply(
        lambda x: x['OrfStarts'].index(int(x['OrfStart'])), axis=1)
    so_df['OrfPosition'] = so_df.apply(lambda x: 'first' if x['OrfIndex'] == 0 else (
        'last' if x['OrfIndex'] == len(x['OrfStarts'])-1 else 'middle'), axis=1)
    return so_df


def count_orfs_by_position(so_df):
    nr_first_orfs = sum(so_df['OrfPosition'] == 'first')
    nr_middle_orfs = sum(so_df['OrfPosition'] == 'middle')
    nr_last_orfs = sum(so_df['OrfPosition'] == 'last')
    return nr_first_orfs, nr_middle_orfs, nr_last_orfs


def val_so_by_position(validated_so_df, nr_validated_so, outdir, region_type):
    validated_so_df = get_so_position_in_transcript(validated_so_df)
    nr_first_orfs, nr_middle_orfs, nr_last_orfs = count_orfs_by_position(
        validated_so_df)

    assert nr_first_orfs + nr_middle_orfs + nr_last_orfs == nr_validated_so

    plot_three_category_pie(nr_first_orfs,
                            nr_middle_orfs,
                            nr_last_orfs,
                            nr_validated_so,
                            ['# first ORFs', '# middle ORFs', '# last ORFs'],
                            'Pie chart of ORF positions of validated Split-ORFs',
                            outdir,
                            'positions_of_validated_SO_pie_chart',
                            region_type,
                            ['#75C1C5', '#FFC500', '#CC79A7']
                            )
    return nr_first_orfs, nr_middle_orfs, nr_last_orfs


def all_URs_by_position(dna_ur_df, all_predicted_so_orfs, outdir, region_type):
    dna_ur_df = pd.merge(
        dna_ur_df, all_predicted_so_orfs.iloc[:, 1:6], on='OrfID', how='left')

    nr_DNA_URs = len(dna_ur_df.index)

    dna_ur_df = get_so_position_in_transcript(dna_ur_df)
    nr_first_orfs, nr_middle_orfs, nr_last_orfs = count_orfs_by_position(
        dna_ur_df)
    plot_three_category_pie(nr_first_orfs,
                            nr_middle_orfs,
                            nr_last_orfs,
                            nr_DNA_URs,
                            ['# first ORFs', '# middle ORFs', '# last ORFs'],
                            'Pie chart of ORF positions of sos with DNA UR',
                            outdir,
                            'positions_of_SO_with_UR_pie_chart',
                            region_type,
                            ['#75C1C5', '#FFC500', '#CC79A7']
                            )
    return nr_first_orfs, nr_middle_orfs, nr_last_orfs


def val_perc_first_middle_last_orfs_csv(nr_val_first_orfs,
                                        nr_val_middle_orfs,
                                        nr_val_last_orfs,
                                        nr_first_orfs_ur,
                                        nr_middle_orfs_ur,
                                        nr_last_orfs_ur,
                                        outdir,
                                        outname):
    perc_val_first_orfs = nr_val_first_orfs/nr_first_orfs_ur
    perc_val_middle_orfs = nr_val_middle_orfs/nr_middle_orfs_ur
    perc_val_last_orfs = nr_val_last_orfs/nr_last_orfs_ur
    perc_val_dict = {'perc_val_first_orfs': [perc_val_first_orfs],
                     'perc_val_middle_orfs': [perc_val_middle_orfs],
                     'perc_val_last_orfs': [perc_val_last_orfs]}
    perc_val_df = pd.DataFrame(perc_val_dict)
    perc_val_df.to_csv(os.path.join(outdir, outname))


def calculate_overlapping_region_percentage(start1, end1, start2, end2):
    if end1 <= start2 or end2 <= start1:
        return 0
    elif start1 < end2 and start2 < end1:
        overlap_start = max(start2, start1)
        overlap_end = min(end1, end2)
        nr_bp_overlap = overlap_end - overlap_start
        shorter_region = min(end2-start2, end1-start1)
        return nr_bp_overlap/shorter_region


def get_max_overlap_of_regions_in_df(chr_df, threshold=0.7):
    for index1 in chr_df.index:
        for index2 in chr_df.index:
            if index1 != index2:
                start1 = float(chr_df.iloc[index1]['start'])
                end1 = float(chr_df.iloc[index1]['stop'])
                start2 = float(chr_df.iloc[index2]['start'])
                end2 = float(chr_df.iloc[index2]['stop'])
                overlap = calculate_overlapping_region_percentage(
                    start1, end1, start2, end2)
                if overlap > threshold:
                    chr_df.iloc[index2]['OrfPositionsOverlapping'].add(
                        chr_df.iloc[index1]['OrfPosition'])
                    chr_df.iloc[index2]['OrfIDsOverlapping'].add(
                        chr_df.iloc[index1]['OrfID'])
                    chr_df.iloc[index1]['OrfPositionsOverlapping'].add(
                        chr_df.iloc[index2]['OrfPosition'])
                    chr_df.iloc[index1]['OrfIDsOverlapping'].add(
                        chr_df.iloc[index2]['OrfID'])
                if overlap > float(chr_df.iloc[index1]['OverlapPercentage']):
                    chr_df.loc[index1, 'OverlapPercentage'] = overlap
                if overlap > float(chr_df.iloc[index2]['OverlapPercentage']):
                    chr_df.loc[index2, 'OverlapPercentage'] = overlap
    return chr_df


def get_transcript_with_2_distinct_val_URs(val_dna_overlapping_ur_df):
    """this information needs to be seen, so printed or saved as a file, otherwise this function is useless"""
    val_dna_overlapping_ur_df['OrfTransID'] = val_dna_overlapping_ur_df['OrfTransID'].apply(
        lambda x: x.split(',')).copy()
    val_dna_overlapping_ur_df['OrfTransID'].explode().value_counts()
    # how often are there two URs that do not overlap more than 50% validated
    sum(val_dna_overlapping_ur_df['OrfTransID'].explode(
    ).value_counts() > 1)
    val_dna_overlapping_ur_df['OrfTransID'].explode().value_counts(
    )[val_dna_overlapping_ur_df['OrfTransID'].explode().value_counts() > 1]


def identify_overlapping_unique_regions(validated_so_df, dna_ur_df, outdir):
    # filter for validated ORFs
    val_dna_ur_df = dna_ur_df[dna_ur_df['OrfID'].isin(
        validated_so_df['OrfID'])].copy()

    # group several exonic URs per ORF together
    val_dna_ur_df = val_dna_ur_df.groupby('OrfID').agg({'start': 'min',
                                                        'stop': 'max',
                                                        'chr': 'first',
                                                        'ID': lambda x: ','.join(x),
                                                        'OrfTransID': 'first'}).reset_index().copy()

    # map ORF positions to IDs
    orf_id_position_map = validated_so_df.set_index('OrfID')['OrfPosition']
    val_dna_ur_df['OrfPosition'] = val_dna_ur_df['OrfID'].map(orf_id_position_map
                                                              )

    # concatenate genomic regions
    val_dna_ur_df['genomic_UR'] = val_dna_ur_df['chr'].astype(
        str) + '_' + val_dna_ur_df['start'].astype(str) + '_' + val_dna_ur_df['stop'].astype(str)

    val_dna_ur_df['OverlapPercentage'] = 0.0
    val_dna_ur_df['OrfPositionsOverlapping'] = val_dna_ur_df['OrfPosition'].apply(
        lambda x: set([x]))
    val_dna_ur_df['OrfIDsOverlapping'] = val_dna_ur_df['OrfID'].apply(
        lambda x: set([x]))

    # get > 70% overlapping URs
    chr_dfs = {chr: chr_df.reset_index(drop=True).copy(
    ) for chr, chr_df in val_dna_ur_df.groupby('chr')}
    chr_dfs = {chr: get_max_overlap_of_regions_in_df(
        chr_df, 0.5) for chr, chr_df in chr_dfs.items()}
    val_dna_overlapping_ur_df = pd.concat(
        chr_dfs.values()).reset_index(drop=True).copy()

    val_dna_overlapping_ur_df['ORFs_sharing_region'] = val_dna_overlapping_ur_df['OrfIDsOverlapping'].apply(
        lambda x: len(x))
    val_dna_overlapping_ur_df['shared_region_type'] = val_dna_overlapping_ur_df['OrfPositionsOverlapping'].apply(
        lambda x: len(x))

    # subset for single ORF having the UR or ORFs from the same position (first, middle last)
    val_dna_overlapping_ur_df = val_dna_overlapping_ur_df[(val_dna_overlapping_ur_df['ORFs_sharing_region'] == 1) | (
        val_dna_overlapping_ur_df['shared_region_type'] == 1)]

    val_dna_overlapping_ur_df.loc[:, 'OrfIDsOverlapping'] = val_dna_overlapping_ur_df['OrfIDsOverlapping'].apply(
        lambda x: frozenset(x))

    val_dna_overlapping_ur_df = val_dna_overlapping_ur_df.groupby('OrfIDsOverlapping').agg(
        {'genomic_UR': 'first',
            'ORFs_sharing_region': 'first',
            'OrfPosition': 'first',
            'ID': lambda x: ','.join(x),
            'OrfTransID': lambda x: ','.join(x),
            'OrfPositionsOverlapping': 'first',
            'OrfIDsOverlapping': 'first',
            'OverlapPercentage': 'max',
         }).reset_index(drop=True)

    get_transcript_with_2_distinct_val_URs(val_dna_overlapping_ur_df)

    val_dna_overlapping_ur_df['OrfPosition'].value_counts(
    ).reset_index().to_csv(os.path.join(outdir, 'distinct_URs_per_position.csv'))

    validated_so_df.to_csv(os.path.join(outdir, 'validated_so_df.csv'))

    return val_dna_overlapping_ur_df
