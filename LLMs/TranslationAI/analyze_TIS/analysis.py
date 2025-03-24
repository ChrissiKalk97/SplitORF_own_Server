import os
import re

import pandas as pd
from scipy import stats
from pybedtools import BedTool


def nr_trans_and_mean_probs(TransAI_SO_preds, TIS_results_df, SO_transcripts):
    print('Number of SO predicted transcripts', len(TransAI_SO_preds.index))
    print('Number of non-SO predicted transcripts',
          len(TIS_results_df.index) - len(TransAI_SO_preds.index))

    SO_2nd_start_mean_prob = TIS_results_df[TIS_results_df['OrfTransID'].isin(
        SO_transcripts)].iloc[:, 2].fillna('0,0').apply(
            lambda x: float(x.split(',')[1])).mean()

    random_2nd_start_mean_prob = TIS_results_df[~TIS_results_df['OrfTransID'].isin(
        SO_transcripts)].iloc[:, 2].fillna('0,0').apply(
            lambda x: float(x.split(',')[1])).mean()

    print('Mean probability of 2nd start in SO transcripts', SO_2nd_start_mean_prob)
    print('Mean probability of 2nd starts in non-SO predicted NMD transcripts',
          random_2nd_start_mean_prob)


def get_TransAI_SO_info(TransAI_SO_preds):
    # How often do positions corresponds?
    TransAI_SO_preds["TIS_dict"] = TransAI_SO_preds[range(1, 11)].apply(
        lambda row: {str(k).split(',')[0]: str(k).split(',')[1] for k in row.dropna()}, axis=1)

    TransAI_SO_preds["TIS_pos_list"] = TransAI_SO_preds[range(1, 11)].apply(
        lambda row: [str(k).split(',')[0] for k in row.dropna()], axis=1)

    return TransAI_SO_preds


def merge_df(TransAI_SO_preds, predicted_SO_ORFs):
    df_merged = pd.merge(TransAI_SO_preds[[
                         "TIS_dict", "TIS_pos_list", "OrfTransID"]], predicted_SO_ORFs, on='OrfTransID', how='outer')

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


def background_probs(TIS_results_df, SO_transcripts):
    background_transcripts_df = TIS_results_df[~TIS_results_df['OrfTransID'].isin(
        SO_transcripts)]
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
    print(f'The best scoring background starts are tested to have a lower distribution than the best scoring SO starts.\n\
          The wilcoxon rank sum test statistic is {result_one_sided[0]}.\n\
            The wilcoxon rank sum p-value is {result_one_sided[1]}')


def get_ORF_start_probs(Ribo_df, TransAI_df, on='OrfTransID'):
    Ribo_df['ORF_start'] = Ribo_df['name'].apply(
        lambda x: int(x.split(':')[2]))
    Ribo_df_merged = pd.merge(
        Ribo_df, TransAI_df, on=on, how='left')
    Ribo_df_merged['ORF_number'] = Ribo_df_merged.apply(
        lambda x: x['OrfStarts_sorted'].index(x['ORF_start']) + 1, axis=1)
    # ['OrfStartsTransAIProbs'] are sorted by the ORF
    Ribo_df_merged['ORF_prob'] = Ribo_df_merged.apply(lambda x: float(x['TIS_dict'][str(
        x['ORF_start'])]) if str(x['ORF_start']) in x['TIS_pos_list'] else 0, axis=1)

    Ribo_df_merged_first = Ribo_df_merged[Ribo_df_merged['ORF_number'] == 1]
    Ribo_df_merged_second = Ribo_df_merged[Ribo_df_merged['ORF_number'] > 1]

    return Ribo_df_merged, Ribo_df_merged_first, Ribo_df_merged_second


def analyze_emp_background_Riboseq(empirical_Ribo_findings_file, df_merged):
    sample = os.path.basename(
        empirical_Ribo_findings_file).rsplit('_', 2)[0]
    print(sample)
    Ribo_results_df = pd.read_csv(
        empirical_Ribo_findings_file, header=0, index_col=0)
    Ribo_results_significant_df = Ribo_results_df[Ribo_results_df["significant"] == 1]
    Ribo_results_significant_df['ORF'] = Ribo_results_significant_df['name'].apply(
        lambda x: x.split(':')[1])
    Ribo_results_significant_df['OrfTransID'] = Ribo_results_significant_df['name'].apply(
        lambda x: x.split(':')[0].split('|')[1])
    Ribo_results_significant_df['genomic_UR'] = Ribo_results_significant_df['new_name'].apply(lambda x: x.split(':')[-1])
    Ribo_results_significant_df = Ribo_results_significant_df.groupby('genomic_UR').agg('first')
    Ribo_results_significant_df = Ribo_results_significant_df.reset_index()


    Ribo_df_merged, Ribo_df_merged_first, Ribo_df_merged_second = get_ORF_start_probs(
        Ribo_results_significant_df, df_merged)
    

    print('Nr validated ORFs not being the first',
          sum(Ribo_df_merged['ORF_number'] > 1))

    return Ribo_df_merged, Ribo_df_merged_first, Ribo_df_merged_second, sample


def create_RiboTISH_BedTool(RiboTISH_SO_df):
    # extract bedformat information
    RiboTISH_SO_df["Chr"] = RiboTISH_SO_df['GenomePos'].apply(
        lambda x: x.split(':')[0])
    RiboTISH_SO_df["Start"] = RiboTISH_SO_df['GenomePos'].apply(
        lambda x: int(re.split(r'[:\-]', x)[1]))
    RiboTISH_SO_df["Stop"] = RiboTISH_SO_df['GenomePos'].apply(
        lambda x: int(re.split(r'[:\-]', x)[2]))
    RiboTISH_SO_df["Strand"] = RiboTISH_SO_df['GenomePos'].apply(
        lambda x: x.split(':')[2])
    RiboTISH_SO_df['Score'] = '0'
    # create a name with all necessary information
    RiboTISH_SO_df["Name"] = RiboTISH_SO_df['Gid'] + '|' + \
        RiboTISH_SO_df['Tid'] + '|' + \
        RiboTISH_SO_df['Symbol'] + '|' + \
        RiboTISH_SO_df['GenomePos'] + '|' + \
        RiboTISH_SO_df['RiboPvalue'].apply(lambda x: str(x)) + '|' + \
        RiboTISH_SO_df['InFrameCount'].apply(lambda x: str(x))
    # create string for BedTools
    BedTool_string = "\n".join(f"{c}\t{s}\t{e}\t{n}\t{sc}\t{st}"
                               for c, s, e, n, sc, st
                               in zip(RiboTISH_SO_df['Chr'],
                                      RiboTISH_SO_df['Start'],
                                      RiboTISH_SO_df['Stop'],
                                      RiboTISH_SO_df['Name'],
                                      RiboTISH_SO_df['Score'],
                                      RiboTISH_SO_df['Strand']))

    RiboTISH_BedTool = BedTool(BedTool_string, from_string=True)
    return RiboTISH_BedTool, RiboTISH_SO_df


def obtain_correct_ORF_RiboTISH(UR_BedTool, RiboTISH_BedTool):
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
        'nr_bp': 'sum'})
    # print('length after aggregation', len(
    #     URs_found_in_RiboTISH_df.index))
    URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df.reset_index()

    # filter for the wrong transcript
    URs_found_in_RiboTISH_df['Tid_RT'] = URs_found_in_RiboTISH_df['rt_name'].apply(
        lambda x: x.split('|')[1])
    URs_found_in_RiboTISH_df['Tid_UR'] = URs_found_in_RiboTISH_df['name'].apply(
        lambda x: re.split(r'[:\|]', x)[1])
    URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df[
        URs_found_in_RiboTISH_df['Tid_RT'] == URs_found_in_RiboTISH_df['Tid_UR']]
    # print('nr regions after  Tid filtering',
    #         len(URs_found_in_RiboTISH_df.index))

    # calculate the UR length
    URs_found_in_RiboTISH_df['UR_length'] = URs_found_in_RiboTISH_df['name'].apply(
        lambda x: int(x.split(':')[5]) - int(x.split(':')[4]))
    # require the overlapping bp to equal the length of the UR
    # even though: is that correct? If I now only have one coordinate min max???
    URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df[
        URs_found_in_RiboTISH_df['UR_length'] == URs_found_in_RiboTISH_df['nr_bp']]
    # print('nr regions after bp overlap filtering',
    #         len(URs_found_in_RiboTISH_df.index))

    idx_max = URs_found_in_RiboTISH_df.groupby('rt_name')[
        'nr_bp'].idxmax()

    # Filter the DataFrame to keep only those rows with the maximum 'col2' values
    URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df.loc[idx_max]
    print('nr regions aftermax bp overlap filtering for dups',
          len(URs_found_in_RiboTISH_df.index))

    return URs_found_in_RiboTISH_df


def check_start_probs_Ensembl_canonical(Ribo_df_merged, Ensembl_canonical_df, verbose = False):
    Ensembl_ribo_merged_df = pd.merge(Ribo_df_merged, Ensembl_canonical_df, on = 'OrfTransID', how = 'left')
    Ensembl_ribo_merged_df['Ensembl_start_in_SO_starts'] = Ensembl_ribo_merged_df.apply(
        lambda x: x['cDNA coding start'] in x['OrfStarts_sorted'], axis = 1)
    Ens_start_not_SO = Ensembl_ribo_merged_df[Ensembl_ribo_merged_df['Ensembl_start_in_SO_starts'] == False]
    Ens_start_not_SO['Ensmebl_start_prob'] = Ens_start_not_SO[['TIS_dict', 'TIS_pos_list', 'cDNA coding start', 'best_start_SO', 'second_best_start_SO']].apply(
        lambda x: x['TIS_dict'][str(x['cDNA coding start'])] if str(x['cDNA coding start']) in x['TIS_pos_list'] else 0, axis = 1)
    
    if verbose:
        print(Ens_start_not_SO[['cDNA coding start',  'OrfStarts', 'Ensmebl_start_prob']])
