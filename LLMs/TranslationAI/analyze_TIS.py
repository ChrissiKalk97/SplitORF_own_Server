# ------------------ IMPORTS ------------------ #
import sys
import os
import re
import glob

import seaborn as sbn
import matplotlib.pyplot as plt
import pandas as pd

from scipy import stats
from pybedtools import BedTool

# ------------------ CONSTANTS ------------------ #
TIS_results = sys.argv[1]
SO_results = sys.argv[2]

outdir = os.path.dirname(TIS_results)
datatype = os.path.basename(SO_results).split('_')[1].split('.')[0]

# TIS_results = '/projects/splitorfs/work/LLMs/TranslationAI/Output/NMD_trnascripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt'
# SO_results = '/projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt'


def main():
    # ------------------ DATA IMPORT ------------------ #
    # only read the first 10 start codons, the rest will not be of interest
    TIS_results_df = pd.read_csv(
        TIS_results, sep="\t", usecols=range(11), header=None)
    TIS_results_df.rename(columns={0: 'fasta_header'}, inplace=True)
    TIS_results_df['OrfTransID'] = TIS_results_df['fasta_header'].apply(
        lambda x: re.split(r'\(|\)|\|', x)[4])

    predicted_SO_ORFs = pd.read_csv(SO_results, header=0, sep='\t')
    SO_transcripts = predicted_SO_ORFs['OrfTransID'].to_list()
    predicted_SO_ORFs['OrfPos'] = predicted_SO_ORFs['OrfPos'].apply(
        lambda x: x.split(','))
    predicted_SO_ORFs['OrfStarts'] = predicted_SO_ORFs['OrfPos'].apply(
        lambda x: [y.split('-')[0] for y in x])
    # get Id to match the TIS transformer output one!
    predicted_SO_ORFs['nr_SO_starts'] = predicted_SO_ORFs['OrfPos'].apply(
        lambda x: len(x))

    # ------------------ COMPARE PREDICTIONS SO AND TRANSLATIONAI ------------------ #
    TransAI_SO_preds = TIS_results_df[TIS_results_df['OrfTransID'].isin(
        SO_transcripts)]
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

    # How often do positions corresponds?
    TransAI_SO_preds["TIS_dict"] = TransAI_SO_preds[range(1, 11)].apply(
        lambda row: {str(k).split(',')[0]: str(k).split(',')[1] for k in row.dropna()}, axis=1)

    TransAI_SO_preds["TIS_pos_list"] = TransAI_SO_preds[range(1, 11)].apply(
        lambda row: [str(k).split(',')[0] for k in row.dropna()], axis=1)

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

    # ------------------ BACKGROUND AUG PROBABILITIES ------------------ #
    # idea: non-SO NMD transcripts as the background, take their 2 best scorign AUGs and
    # compare to the 2 best scroing SO
    # compare the best with the best and the second best with the second best

    # the first column will have the best scoring starts and the second column the second best scoring start
    background_transcripts_df = TIS_results_df[~TIS_results_df['OrfTransID'].isin(
        SO_transcripts)]
    background_transcripts_df['best_start'] = background_transcripts_df[1].fillna(
        '0,0').apply(lambda x: float(x.split(',')[1]))
    background_transcripts_df['second_best_start'] = background_transcripts_df[2].fillna(
        '0,0').apply(lambda x: float(x.split(',')[1]))

    # with the following the information of which start belongs to which probability
    # will be lost, hence make a new column and keep the old one
    df_merged['OrfStartsTransAIProbs_sorted'] = df_merged['OrfStartsTransAIProbs'].apply(
        lambda x: sorted([float(start) for start in x], reverse=True))
    df_merged['best_start_SO'] = df_merged['OrfStartsTransAIProbs_sorted'].apply(
        lambda x: x[0])
    df_merged['second_best_start_SO'] = df_merged['OrfStartsTransAIProbs_sorted'].apply(
        lambda x: x[1])

    # ------------------ PLOT BEST PROBABILITIES ------------------ #
    plt.figure(figsize=(8, 6))
    sbn.histplot(df_merged['best_start_SO'], color="blue",
                 label="SO best start", bins=10, alpha=0.6)
    sbn.histplot(background_transcripts_df['best_start'],
                 color="red", label="Non-SO best start", bins=10, alpha=0.6)

    # Labels and legend
    plt.xlabel("Probability of the best start site")
    plt.ylabel("Frequency")
    plt.title(
        f"Distribution of Best Start Probabilities for SO and non-SO {datatype} transcripts")
    plt.legend()

    # Save plot
    plt.savefig(f"{outdir}/{datatype}_best_start_SO_non_SO_historgam.svg",
                dpi=300, bbox_inches="tight")
    plt.close()

    # ------------------ WILOCXON RANK SUM TEST FOR THE BEST STARTS ------------------ #
    # one sided: interested whether the random starts are less likely than the SO starts
    result_one_sided = stats.ranksums(
        background_transcripts_df['best_start'], df_merged['best_start_SO'], alternative='less')
    print(f'The best scoring background starts are tested to have a lower distribution than the best scoring SO starts.\n\
          The wilcoxon rank sum test statistic is {result_one_sided[0]}.\n\
            The wilcoxon rank sum p-value is {result_one_sided[1]}')
    # stats.ranksums(background_transcripts_df['best_start'], df_merged['best_start_SO'], alternative='two-sided')

    # ------------------ PLOT SECOND BEST PROBABILITIES ------------------ #
    plt.figure(figsize=(8, 6))
    sbn.histplot(df_merged['second_best_start_SO'], color="blue",
                 label="SO second best start", bins=10, alpha=0.6)
    sbn.histplot(background_transcripts_df['second_best_start'],
                 color="red", label="Non-SO second best start", bins=10, alpha=0.6)

    # Labels and legend
    plt.xlabel("Probability of the second best start site")
    plt.ylabel("Frequency")
    plt.title(
        f"Distribution of Second Best Start Probabilities for SO and non-SO {datatype} transcripts")
    plt.legend()

    # Save plot
    plt.savefig(f"{outdir}/{datatype}_second_best_start_SO_non_SO_historgam.svg",
                dpi=300, bbox_inches="tight")
    plt.close()

    # ------------------ WILOCXON RANK SUM TEST FOR THE BEST STARTS ------------------ #
    # one sided: interested whether the random starts are less likely than the SO starts
    result_one_sided = stats.ranksums(
        background_transcripts_df['second_best_start'], df_merged['second_best_start_SO'], alternative='less')
    print(f'The best scoring background starts are tested to have a lower distribution than the best scoring SO starts.\n\
          The wilcoxon rank sum test statistic is {result_one_sided[0]}.\n\
            The wilcoxon rank sum p-value is {result_one_sided[1]}')

    #################################################################################
    # ------------------ COMPARE WITH RIBOSEQ EMPIRICAL FINDINGS ------------------ #
    #################################################################################
    df_merged['OrfStarts'] = df_merged['OrfStarts'].apply(
        lambda x: [int(start) for start in x])
    df_merged['OrfStarts_sorted'] = df_merged['OrfStarts'].apply(
        lambda x: sorted(x))

    for empirical_Ribo_findings_file in glob.glob("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample/NMD_genome/*_unique_regions.csv"):
        sample = os.path.basename(
            empirical_Ribo_findings_file).rsplit('_', 2)[0]
        print(sample)
        Ribo_results_df = pd.read_csv(
            empirical_Ribo_findings_file, header=0, index_col=0)
        Ribo_results_significant_df = Ribo_results_df[Ribo_results_df["significant"] == 1]
        Ribo_results_significant_df['ORF'] = Ribo_results_df['name'].apply(
            lambda x: x.split(':')[1])
        Ribo_results_significant_df['ORF_start'] = Ribo_results_df['name'].apply(
            lambda x: int(x.split(':')[2]))
        Ribo_results_significant_df['OrfTransID'] = Ribo_results_df['name'].apply(
            lambda x: x.split(':')[0].split('|')[1])
        Ribo_df_merged = pd.merge(
            Ribo_results_significant_df, df_merged, on='OrfTransID', how='left')
        Ribo_df_merged['ORF_number'] = Ribo_df_merged.apply(
            lambda x: x['OrfStarts_sorted'].index(x['ORF_start']) + 1, axis=1)
        # ['OrfStartsTransAIProbs'] are sorted by the ORF
        Ribo_df_merged['ORF_prob'] = Ribo_df_merged.apply(lambda x: float(x['TIS_dict'][str(
            x['ORF_start'])]) if str(x['ORF_start']) in x['TIS_pos_list'] else 0, axis=1)

        Ribo_df_merged_first = Ribo_df_merged[Ribo_df_merged['ORF_number'] == 1]
        Ribo_df_merged_second = Ribo_df_merged[Ribo_df_merged['ORF_number'] > 1]

        print(sum(Ribo_df_merged['ORF_number'] > 1))
        # 57 from 300 thats a lot!!!

        plt.figure(figsize=(8, 6))
        sbn.histplot(Ribo_df_merged_first['ORF_prob'], color="blue",
                     label="First val ORFs", bins=10, alpha=0.6)
        sbn.histplot(Ribo_df_merged_second['ORF_prob'], color="red",
                     label="Second val ORFs", bins=10, alpha=0.6)

        # Labels and legend
        plt.xlabel("Probability of first and second Riboseq val ORFs")
        plt.ylabel("Frequency")
        plt.title(
            f"Distribution of first and second Riboseq val ORF probs - {datatype}")
        plt.legend()

        plt.savefig(f"{outdir}/{datatype}_{sample}_Riboseq_empirical_first_and_second_ORF.svg",
                    dpi=300, bbox_inches="tight")
        plt.close()

    #################################################################################
    # ------------------ COMPARE WITH RIBOTISH VAL ORFS FINDINGS ------------------ #
    #################################################################################
    # create unique genomic region bedtool
    UR_BedTool = BedTool(
        f'/projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_{datatype}_16_12_24_chrom_sorted.bed')
    
    # get TranslationAI information required
    df_merged_RiboTISH = df_merged[['TIS_dict', 'OrfTransID', 'TIS_pos_list', 'OrfStarts_sorted']]
    df_merged_RiboTISH = df_merged_RiboTISH.rename(columns = {'OrfTransID': 'Tid_RT'})


    # Reformat to bed format, pybedtools
    for RiboTISH_file in glob.glob("/projects/splitorfs/work/Riboseq/Output/RiboTISH_NMD_custom/OHMX20220060_001*.csv"):
        RiboTISH_SO_df = pd.read_csv(RiboTISH_file, header=0)
        # extract bedformat information
        RiboTISH_SO_df["Chr"] = RiboTISH_SO_df['GenomePos'].apply(
            lambda x: int(x.split(':')[0]))
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
        print('length before aggregation', len(URs_found_in_RiboTISH_df.index))
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
        print('length after aggregation', len(URs_found_in_RiboTISH_df.index))
        URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df.reset_index()

        # filter for the wrong transcript
        URs_found_in_RiboTISH_df['Tid_RT'] = URs_found_in_RiboTISH_df['rt_name'].apply(lambda x: x.split('|')[1])
        URs_found_in_RiboTISH_df['Tid_UR'] = URs_found_in_RiboTISH_df['name'].apply(lambda x: re.split(r'[:\|]', x)[1])
        URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df[URs_found_in_RiboTISH_df['Tid_RT'] == URs_found_in_RiboTISH_df['Tid_UR']]
        print('nr regions after  Tid filtering', len(URs_found_in_RiboTISH_df.index))

        # calculate the UR length
        URs_found_in_RiboTISH_df['UR_length'] = URs_found_in_RiboTISH_df['name'].apply(
            lambda x: int(x.split(':')[5]) - int(x.split(':')[4]))
        # require the overlapping bp to equal the length of the UR
        # even though: is that correct? If I now only have one coordinate min max???
        URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df[
            URs_found_in_RiboTISH_df['UR_length'] == URs_found_in_RiboTISH_df['nr_bp']]
        print('nr regions after bp overlap filtering', len(URs_found_in_RiboTISH_df.index))

        idx_max = URs_found_in_RiboTISH_df.groupby('rt_name')['nr_bp'].idxmax()

        # Filter the DataFrame to keep only those rows with the maximum 'col2' values
        URs_found_in_RiboTISH_df = URs_found_in_RiboTISH_df.loc[idx_max]
        print('nr regions aftermax bp overlap filtering for dups', len(URs_found_in_RiboTISH_df.index))

        # Get transcript coordinates ORF starts
        URs_found_in_RiboTISH_df['ORF_start'] = URs_found_in_RiboTISH_df['name'].apply(lambda x: int(x.split(':')[2]))

        URs_RiboTISH_TransAI_df = pd.merge(URs_found_in_RiboTISH_df, df_merged_RiboTISH, on = 'Tid_RT', how = 'left')
        URs_RiboTISH_TransAI_df['ORF_prob'] = URs_RiboTISH_TransAI_df.apply(lambda x: float(x['TIS_dict'][str(
            x['ORF_start'])]) if str(x['ORF_start']) in x['TIS_pos_list'] else 0, axis=1)
        URs_RiboTISH_TransAI_df['ORF_number'] = URs_RiboTISH_TransAI_df.apply(
            lambda x: x['OrfStarts_sorted'].index(x['ORF_start']) + 1, axis=1)
        
        URs_RiboTISH_TransAI_df_first = URs_RiboTISH_TransAI_df[URs_RiboTISH_TransAI_df['ORF_number'] == 1]
        URs_RiboTISH_TransAI_df_second = URs_RiboTISH_TransAI_df[URs_RiboTISH_TransAI_df['ORF_number'] > 1]
        
        plt.figure(figsize=(8, 6))
        sbn.histplot(URs_RiboTISH_TransAI_df_first['ORF_prob'], color="blue",
                     label="First val ORFs", bins=10, alpha=0.6)
        sbn.histplot(URs_RiboTISH_TransAI_df_second['ORF_prob'], color="red",
                     label="Second val ORFs", bins=10, alpha=0.6)

        # Labels and legend
        plt.xlabel("Probability of first and second RiboTISH val ORFs")
        plt.ylabel("Frequency")
        plt.title(
            f"Distribution of first and second RiboTISH val ORF probs - {datatype}")
        plt.legend()

        plt.show()



if __name__ == "__main__":
    main()
