import os
import re

import pandas as pd
from scipy import stats
from pybedtools import BedTool

from plotting import plot_three_category_pie, plot_so_start_probs_by_position_box

#################################################################################
# ------------------ ALL PRED SO AGAISNT BACKGROUND          ------------------ #
#################################################################################


def nr_trans_and_mean_probs(TransAI_so_preds, TIS_results_df, so_transcripts):
    print('Number of SO predicted transcripts', len(TransAI_so_preds.index))
    print('Number of non-SO predicted transcripts',
          len(TIS_results_df.index) - len(TransAI_so_preds.index))

    so_2nd_start_mean_prob = TIS_results_df[TIS_results_df['OrfTransID'].isin(
        so_transcripts)].iloc[:, 2].fillna('0,0').apply(
            lambda x: float(x.split(',')[1])).mean()

    random_2nd_start_mean_prob = TIS_results_df[~TIS_results_df['OrfTransID'].isin(
        so_transcripts)].iloc[:, 2].fillna('0,0').apply(
            lambda x: float(x.split(',')[1])).mean()

    print('Mean probability of 2nd start in SO transcripts', so_2nd_start_mean_prob)
    print('Mean probability of 2nd starts in non-SO predicted NMD transcripts',
          random_2nd_start_mean_prob)


def get_TransAI_so_info(TransAI_so_preds):
    # How often do positions corresponds?
    TransAI_so_preds["TIS_dict"] = TransAI_so_preds[range(1, 11)].apply(
        lambda row: {str(k).split(',')[0]: str(k).split(',')[1] for k in row.dropna()}, axis=1)

    TransAI_so_preds["TIS_pos_list"] = TransAI_so_preds[range(1, 11)].apply(
        lambda row: [str(k).split(',')[0] for k in row.dropna()], axis=1)

    return TransAI_so_preds


def merge_df(TransAI_so_preds, predicted_so_orfs):
    df_merged = pd.merge(TransAI_so_preds[[
                         "TIS_dict", "TIS_pos_list", "OrfTransID"]], predicted_so_orfs, on='OrfTransID', how='outer')

    # Are k ORFs predicted in the top k?
    # First sample the TranslationAI predicted positions for the same number as the predicted SOs for that transcript
    df_merged['top_k_TIS'] = df_merged.apply(
        lambda x: x['TIS_pos_list'][:len(x['OrfStarts'])], axis=1)
    df_merged['SOs_are_top_TIS'] = df_merged.apply(lambda x: all(
        start in x['top_k_TIS'] for start in x['OrfStarts']), axis=1)

    print('Number of times that all k SO starts are within the topk predicted TIS by TranslationaAI',
          df_merged['SOs_are_top_TIS'].sum())

    # Find out how often at least 2 predicted SO start positions obtain a probability of higher than or equal to 10%
    # this will get the predicted probabilities for the SO start sites in the order that they are listed in 'OrfStarts'
    df_merged['OrfStartsTransAIProbs'] = df_merged.apply(
        lambda x: [x['TIS_dict'][start] if start in x['TIS_pos_list'] else 0 for start in x['OrfStarts']], axis=1)
    df_merged['NrOrfStartsGreater0_1'] = df_merged['OrfStartsTransAIProbs'].apply(
        lambda x: len([prob for prob in x if float(prob) >= 0.1]))
    print("Number of times at least 2 start codons of the SO pipeline have at least 0.1 prediction probability by TranslationAI: ", sum(
        df_merged['NrOrfStartsGreater0_1'] > 1))

    # with the following the information of which start belongs to which probability
    # will be lost, hence make a new column and keep the old one
    df_merged['OrfStartsTransAIProbs_sorted'] = df_merged['OrfStartsTransAIProbs'].apply(
        lambda x: sorted([float(start) for start in x], reverse=True))
    df_merged['best_start_SO'] = df_merged['OrfStartsTransAIProbs_sorted'].apply(
        lambda x: x[0])
    df_merged['second_best_start_SO'] = df_merged['OrfStartsTransAIProbs_sorted'].apply(
        lambda x: x[1])
    df_merged['OrfStarts'] = df_merged['OrfStarts'].apply(
        lambda x: [int(start) for start in x])
    df_merged['OrfStarts_sorted'] = df_merged['OrfStarts'].apply(
        lambda x: sorted(x))
    return df_merged


def background_probs(TIS_results_df, so_transcripts):
    background_transcripts_df = TIS_results_df[~TIS_results_df['OrfTransID'].isin(
        so_transcripts)].copy()
    # columns are called 1-10 and this is the order of their pred probabilities
    background_transcripts_df['best_start'] = background_transcripts_df[1].fillna(
        '0,0').apply(lambda x: float(x.split(',')[1]))
    background_transcripts_df['second_best_start'] = background_transcripts_df[2].fillna(
        '0,0').apply(lambda x: float(x.split(',')[1]))
    return background_transcripts_df


def wilcoxon_rank_sums(background_transcripts_df, df_merged, nr_start):
    # one sided: interested whether the random starts are less likely than the SO starts
    result_one_sided = stats.ranksums(
        background_transcripts_df[f'{nr_start}_start'], df_merged[f'{nr_start}_start_SO'], alternative='less')
    print(f'The {nr_start} scoring background starts are tested to have a lower distribution than the {nr_start} scoring SO starts.\n\
          The wilcoxon rank sum test statistic is {result_one_sided[0]}.\n\
            The wilcoxon rank sum p-value is {result_one_sided[1]}')


def get_orf_start_probs(Ribo_df, TransAI_df, on='OrfTransID'):
    if 'OrfStart' not in Ribo_df.columns:
        Ribo_df['OrfStart'] = Ribo_df['name'].apply(
            lambda x: int(x.split(':')[2]))
    Ribo_df_merged = pd.merge(
        Ribo_df, TransAI_df, on=on, how='left')
    if 'OrfNumber' not in Ribo_df.columns:
        Ribo_df_merged['OrfNumber'] = Ribo_df_merged.apply(
            lambda x: x['OrfStarts_sorted'].index(x['OrfStart']) + 1, axis=1)
    # ['OrfStartsTransAIProbs'] are sorted by the ORF
    Ribo_df_merged['OrfProb'] = Ribo_df_merged.apply(lambda x: float(x['TIS_dict'][str(
        x['OrfStart'])]) if str(x['OrfStart']) in x['TIS_pos_list'] else 0, axis=1)

    Ribo_df_merged_first = Ribo_df_merged[Ribo_df_merged['OrfNumber'] == 1].copy(
    )
    Ribo_df_merged_second = Ribo_df_merged[Ribo_df_merged['OrfNumber'] > 1].copy(
    )

    return Ribo_df_merged, Ribo_df_merged_first, Ribo_df_merged_second


#################################################################################
# ------------------ FUNCTIONS SO VALIDATED ANALYSIS         ------------------ #
#################################################################################

def analyze_emp_background_riboseq(empirical_Ribo_findings_file, df_merged):
    sample = os.path.basename(
        empirical_Ribo_findings_file).rsplit('_', 2)[0]
    Ribo_results_df = pd.read_csv(
        empirical_Ribo_findings_file, header=0, index_col=0)
    Ribo_results_significant_df = Ribo_results_df[Ribo_results_df["significant"] == 1].copy(
    )
    Ribo_results_significant_df['ORF'] = Ribo_results_significant_df['name'].apply(
        lambda x: x.split(':')[1])
    Ribo_results_significant_df['OrfTransID'] = Ribo_results_significant_df['name'].apply(
        lambda x: x.split(':')[0].split('|')[1])
    Ribo_results_significant_df['genomic_UR'] = Ribo_results_significant_df['new_name'].apply(
        lambda x: x.split(':')[-1])
    # in case an ORF has several URs, select the first entry
    Ribo_results_significant_df = Ribo_results_significant_df.groupby(
        'ORF').agg('first').copy()
    Ribo_results_significant_df = Ribo_results_significant_df.reset_index()
    Ribo_df_merged, Ribo_df_merged_first, Ribo_df_merged_second = get_orf_start_probs(
        Ribo_results_significant_df, df_merged)

    return Ribo_df_merged, Ribo_df_merged_first, Ribo_df_merged_second, sample


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


def print_validation_percentage(nr_val, nr_total):
    val_perc = nr_val/nr_total
    print(f'{val_perc :.2f}')


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


def explode_so_df(predicted_so_orfs):
    predicted_so_orfs = predicted_so_orfs[[
        'OrfTransID', 'OrfIDs', 'OrfStarts']].copy()
    predicted_so_orfs['OrfID'] = predicted_so_orfs.apply(
        lambda x: x['OrfIDs'].split(','), axis=1)
    predicted_so_orfs['OrfStart'] = predicted_so_orfs['OrfStarts']
    all_predicted_so_orfs = predicted_so_orfs.explode(
        ['OrfID', 'OrfStart'], ignore_index=True).copy()

    assert len(all_predicted_so_orfs['OrfTransID'].unique()) == 9261 or len(
        all_predicted_so_orfs['OrfTransID'].unique()) == 5625

    total_nr_so = len(all_predicted_so_orfs.index)
    total_nr_so_transcripts = len(all_predicted_so_orfs['OrfTransID'].unique())
    print('Total number of predicted Split-ORFs', total_nr_so)
    print('Total number of predicted Split-ORF transcripts',
          total_nr_so_transcripts)

    return all_predicted_so_orfs, predicted_so_orfs, total_nr_so


def subset_validated_sos_df(all_predicted_so_orfs):
    all_predicted_so_orfs['ValidationCount'] = all_predicted_so_orfs.iloc[:, range(
        5, 23)].sum(axis=1, numeric_only=True)
    validated_so_df = all_predicted_so_orfs[all_predicted_so_orfs['ValidationCount'] > 1].copy(
    )

    nr_validated_so = len(validated_so_df.index)
    nr_validated_transcripts = len(validated_so_df['OrfTransID'].unique())

    print('Number of validated SOs:', nr_validated_so)
    print('Number of validated SO transcripts:', nr_validated_transcripts)

    return validated_so_df, nr_validated_so, nr_validated_transcripts


def subset_UR_for_expressed_genes(genes_to_keep_file, dna_ur_df, validated_so_df, df_merged):
    # load genes to keep
    genes_to_keep = pd.read_csv(genes_to_keep_file, header=None)

    # get gene and transcript ID in dna_ur_df from name
    dna_ur_df['gene'] = dna_ur_df['OrfTransID'].apply(
        lambda x: x.split('|')[0])
    dna_ur_df['tID'] = dna_ur_df['OrfTransID'].apply(
        lambda x: x.split('|')[1])

    # subset and count ORFs with URs for expressed genes
    dna_ur_df_subset = dna_ur_df[(dna_ur_df['gene'].isin(genes_to_keep[0])) | (
        dna_ur_df['tID'].isin(validated_so_df['OrfTransID']))]
    nr_orfs_with_UR_expressed = len(dna_ur_df_subset['OrfID'].unique())
    # nr_orfs_with_UR = len(dna_ur_df['OrfID'].unique())

    # test whether none of the validated SOs would be lost
    # validated_so_with_gene = df_merged[df_merged['OrfTransID'].isin(
    #     validated_so_df['OrfTransID'])]
    # nr_val_SO_trans = len(validated_so_with_gene.index)
    # nr_val_SO_trans_filtered = len(
    #     validated_so_with_gene[validated_so_with_gene['geneID'].isin(genes_to_keep[0])].index)
    # nr_val_lost = nr_val_SO_trans - nr_val_SO_trans_filtered
    return nr_orfs_with_UR_expressed, dna_ur_df_subset


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
                            ['#75C1C5', '#CC79A7', '#FFC500']
                            )
    return nr_first_orfs, nr_middle_orfs, nr_last_orfs


def all_URs_by_position(dna_ur_df, all_predicted_so_orfs, outdir, region_type):
    dna_ur_df = pd.merge(
        dna_ur_df, all_predicted_so_orfs.iloc[:, 1:5], on='OrfID', how='left')

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
                            'positions_of_so_with_UR_pie_chart',
                            region_type,
                            ['#75C1C5', '#CC79A7', '#FFC500']
                            )
    return nr_first_orfs, nr_middle_orfs, nr_last_orfs


def get_trans_ai_so_preds_df(TransAI_so_preds, validated_so_df):
    TransAI_so_preds_subset = TransAI_so_preds.loc[:, [
        'OrfTransID', 'TIS_dict', 'TIS_pos_list']].copy()
    validated_so_df_subset = validated_so_df.loc[:, [
        'OrfTransID', 'OrfIDs', 'OrfStarts', 'OrfStart', 'OrfID', 'OrfIndex', 'OrfPosition']].copy()
    validated_so_df_subset['OrfNumber'] = validated_so_df_subset['OrfIndex'] + 1

    val_so_trans_ai_df, _, _ = get_orf_start_probs(
        validated_so_df_subset, TransAI_so_preds_subset, on='OrfTransID')
    return val_so_trans_ai_df


# def subset_so_df_by_position(so_df):
#     first_orfs_df = so_df[so_df['OrfPosition'] == 'first'].copy()
#     middle_orfs_df = so_df[so_df['OrfPosition'] == 'middle'].copy()
#     last_orfs_df = so_df[so_df['OrfPosition'] == 'last'].copy()
#     return first_orfs_df, middle_orfs_df, last_orfs_df

# def test_differences_in_start_probs_by_so_position_and_class(all_so_trans_ai_df):
#     first_orfs_df, middle_orfs_df, last_orfs_df = subset_so_df_by_position(all_so_trans_ai_df)

#     first_orfs_df['SO_class']
#     result_one_sided = stats.ranksums(
#     background_transcripts_df[f'{nr_start}_start'], df_merged[f'{nr_start}_start_SO'], alternative='less')


def compare_so_set_probabilities_by_position(all_so_trans_ai_df, dna_ur_df, val_so_trans_ai_df, region_type, outdir, gene_filtering=False):
    if gene_filtering:
        all_so_trans_ai_df['SO_class'] = 'SO without UR or gene not expressed'
    else:
        all_so_trans_ai_df['SO_class'] = 'SO without UR'
    so_with_ur_list = dna_ur_df['OrfID'].to_list()
    so_validated_list = val_so_trans_ai_df['OrfID'].to_list()
    so_with_ur_not_validated_list = [
        so for so in so_with_ur_list if so not in so_validated_list]

    all_so_trans_ai_df.loc[all_so_trans_ai_df['OrfID'].isin(
        so_validated_list), 'SO_class'] = 'validated'
    if gene_filtering:
        all_so_trans_ai_df.loc[all_so_trans_ai_df['OrfID'].isin(
            so_with_ur_not_validated_list), 'SO_class'] = 'SO with UR not validated but gene expressed'
    else:
        all_so_trans_ai_df.loc[all_so_trans_ai_df['OrfID'].isin(
            so_with_ur_not_validated_list), 'SO_class'] = 'SO with UR not validated'

    plot_so_start_probs_by_position_box(
        all_so_trans_ai_df, 'Boxplot of start probabilities by SO position and SO class', region_type, outdir)


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
    sum(val_dna_overlapping_ur_df['OrfTransID'].explode().value_counts() > 1)
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

    return val_dna_overlapping_ur_df

#################################################################################
# ------------------ ENSMEBL CANONICAL                       ------------------ #
#################################################################################


def check_start_probs_Ensembl_canonical(Ribo_df_merged, Ensembl_canonical_df, verbose=False):
    Ensembl_ribo_merged_df = pd.merge(
        Ribo_df_merged, Ensembl_canonical_df, on='OrfTransID', how='left')
    Ensembl_ribo_merged_df['Ensembl_start_in_SO_starts'] = Ensembl_ribo_merged_df.apply(
        lambda x: x['cDNA coding start'] in x['OrfStarts_sorted'], axis=1)
    Ens_start_not_so = Ensembl_ribo_merged_df[Ensembl_ribo_merged_df['Ensembl_start_in_SO_starts'] == False].copy(
    )
    Ens_start_not_so['Ensmebl_start_prob'] = Ens_start_not_so[['TIS_dict', 'TIS_pos_list', 'cDNA coding start', 'best_start_SO', 'second_best_start_SO']].apply(
        lambda x: x['TIS_dict'][str(x['cDNA coding start'])] if str(x['cDNA coding start']) in x['TIS_pos_list'] else 0, axis=1)

    if verbose:
        print(Ens_start_not_so[['cDNA coding start',
              'OrfStarts', 'Ensmebl_start_prob']])


def ensembl_nmd_starts_are_so_starts(Ensembl_canonical_df, all_so_trans_ai_df, val_so_trans_ai_df):
    pass


#################################################################################
# ------------------ RIBOTISH FUNCTIONS                      ------------------ #
#################################################################################

def preprocess_RiboTISH_files(df_merged, datatype, UR_path):
    # create unique genomic region bedtool
    UR_BedTool = BedTool(
        UR_path)

    # get TranslationAI information required
    df_merged_RiboTISH = df_merged[[
        'TIS_dict', 'OrfTransID', 'TIS_pos_list', 'OrfStarts_sorted']].copy()
    df_merged_RiboTISH = df_merged_RiboTISH.rename(
        columns={'OrfTransID': 'Tid_RT'})

    return UR_BedTool, df_merged_RiboTISH


def create_RiboTISH_BedTool(RiboTISH_so_df):
    # extract bedformat information
    RiboTISH_so_df["Chr"] = RiboTISH_so_df['GenomePos'].apply(
        lambda x: x.split(':')[0])
    RiboTISH_so_df["Start"] = RiboTISH_so_df['GenomePos'].apply(
        lambda x: int(re.split(r'[:\-]', x)[1]))
    RiboTISH_so_df["Stop"] = RiboTISH_so_df['GenomePos'].apply(
        lambda x: int(re.split(r'[:\-]', x)[2]))
    RiboTISH_so_df["Strand"] = RiboTISH_so_df['GenomePos'].apply(
        lambda x: x.split(':')[2])
    RiboTISH_so_df['Score'] = '0'
    # create a name with all necessary information
    RiboTISH_so_df["Name"] = RiboTISH_so_df['Gid'] + '|' + \
        RiboTISH_so_df['Tid'] + '|' + \
        RiboTISH_so_df['Symbol'] + '|' + \
        RiboTISH_so_df['GenomePos'] + '|' + \
        RiboTISH_so_df['RiboPvalue'].apply(lambda x: str(x)) + '|' + \
        RiboTISH_so_df['InFrameCount'].apply(lambda x: str(x))
    # create string for BedTools
    BedTool_string = "\n".join(f"{c}\t{s}\t{e}\t{n}\t{sc}\t{st}"
                               for c, s, e, n, sc, st
                               in zip(RiboTISH_so_df['Chr'],
                                      RiboTISH_so_df['Start'],
                                      RiboTISH_so_df['Stop'],
                                      RiboTISH_so_df['Name'],
                                      RiboTISH_so_df['Score'],
                                      RiboTISH_so_df['Strand']))

    RiboTISH_BedTool = BedTool(BedTool_string, from_string=True)
    return RiboTISH_BedTool, RiboTISH_so_df


def obtain_correct_orf_RiboTISH(UR_BedTool, RiboTISH_BedTool):
    # unique region need to be fully covered by modified URs (RiboTISH, f=1.0)
    URs_found_in_RiboTISH_df = UR_BedTool.intersect(
        RiboTISH_BedTool, s=True, f=1.0, wo=True).to_dataframe(disable_auto_names=True, header=None)
    URs_found_in_RiboTISH_df.columns = ['chr',
                                        'start',
                                        'stop',
                                        'name',
                                        'score',
                                        'strand',
                                        'rt_chr',
                                        'rt_start',
                                        'rt_stop',
                                        'rt_name',
                                        'rt_score',
                                        'rt_strand',
                                        'nr_bp']
    # print('length before aggregation', len(
    #     URs_found_in_RiboTISH_df.index))
    URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df.groupby(['name', 'rt_name']).agg({
        'chr': 'first',
        'start': 'min',
        'stop': 'max',
        'strand': 'first',
        'rt_chr': 'first',
        'rt_start': 'first',
        'rt_stop': 'first',
        'rt_strand': 'first',
        'nr_bp': 'sum'}).copy()

    URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df.reset_index()

    # filter for the wrong transcript
    URs_found_in_RiboTISH_df['Tid_RT'] = URs_found_in_RiboTISH_df['rt_name'].apply(
        lambda x: x.split('|')[1])
    URs_found_in_RiboTISH_df['Tid_UR'] = URs_found_in_RiboTISH_df['name'].apply(
        lambda x: re.split(r'[:\|]', x)[1])
    URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df[
        URs_found_in_RiboTISH_df['Tid_RT'] == URs_found_in_RiboTISH_df['Tid_UR']]

    # calculate the UR length
    URs_found_in_RiboTISH_df['UR_length'] = URs_found_in_RiboTISH_df['name'].apply(
        lambda x: int(x.split(':')[5]) - int(x.split(':')[4]))
    # require the overlapping bp to equal the length of the UR
    # even though: is that correct? If I now only have one coordinate min max???
    URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df[
        URs_found_in_RiboTISH_df['UR_length'] == URs_found_in_RiboTISH_df['nr_bp']].copy()

    idx_max = URs_found_in_RiboTISH_df.groupby('rt_name')[
        'nr_bp'].idxmax()

    # Filter the DataFrame to keep only those rows with the maximum 'col2' values
    URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df.loc[idx_max].copy()

    return URs_found_in_RiboTISH_df


#################################################################################
# ------------------ k4neo FUNCTIONS                         ------------------ #
#################################################################################
def preprocess_k4neo_data(k4neo_path, UR_BedTool, df_merged):
    k4neo_validated_transcripts = pd.read_csv(
        k4neo_path, index_col=None, header=None, names=['OrfTransID'])
    k4neo_validated_transcripts = k4neo_validated_transcripts['OrfTransID'].to_list(
    )

    # get URs
    UR_names = []
    for interval in UR_BedTool:
        UR_names.append(interval.name)

    UR_names = [re.split(r'[:|]', name)[1] for name in UR_names]

    # filter k4neo transcripts for transcripts that still have a UR predicted with the
    # current pipeline
    k4neo_validated_transcripts = [
        trans for trans in k4neo_validated_transcripts if trans in UR_names]

    k4neo_TransAI_preds = df_merged[df_merged['OrfTransID'].isin(
        k4neo_validated_transcripts)].copy()

    k4neo_TransAI_preds['best_start'] = k4neo_TransAI_preds['best_start_SO']
    k4neo_TransAI_preds['second_best_start'] = k4neo_TransAI_preds['second_best_start_SO']

    return k4neo_TransAI_preds, k4neo_validated_transcripts
