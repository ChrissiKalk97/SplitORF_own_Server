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

from data_loader import load_TranslationAI, load_SO_results
from analysis import nr_trans_and_mean_probs, get_TransAI_SO_info, merge_df, \
    background_probs, wilcoxon_rank_sums, analyze_emp_background_Riboseq, \
    create_RiboTISH_BedTool, obtain_correct_ORF_RiboTISH, get_ORF_start_probs
from plotting import plot_SO_background, plot_emp_background_TransAI, plot_RiboTISH_TransAI

# ------------------ CONSTANTS ------------------ #
TIS_results = sys.argv[1]
SO_results = sys.argv[2]

outdir = os.path.dirname(TIS_results)
datatype = os.path.basename(SO_results).split('_')[1].split('.')[0]

# TIS_results = '/projects/splitorfs/work/LLMs/TranslationAI/Output/NMD_trnascripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt'
# SO_results = '/projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt'


def main():
    # ------------------ DATA IMPORT ------------------ #
    TIS_results_df = load_TranslationAI(TIS_results)
    predicted_SO_ORFs, SO_transcripts = load_SO_results(SO_results)

    # ------------------ COMPARE PREDICTIONS SO AND TRANSLATIONAI ------------------ #
    TransAI_SO_preds = TIS_results_df[TIS_results_df['OrfTransID'].isin(
        SO_transcripts)]
    nr_trans_and_mean_probs(TransAI_SO_preds, TIS_results_df, SO_transcripts)

    TransAI_SO_preds = get_TransAI_SO_info(TransAI_SO_preds)

    df_merged = merge_df(TransAI_SO_preds, predicted_SO_ORFs)

    # ------------------ BACKGROUND AUG PROBABILITIES ------------------ #
    # idea: non-SO NMD transcripts as the background, take their 2 best scorign AUGs and
    # compare to the 2 best scroing SO
    # compare the best with the best and the second best with the second best

    background_transcripts_df = background_probs(
        TIS_results_df, SO_transcripts)

    # ------------------ PLOT BEST PROBABILITIES ------------------ #
    plot_SO_background(df_merged, background_transcripts_df,
                       'best', datatype, outdir)

    # ------------------ WILOCXON RANK SUM TEST FOR THE BEST STARTS ------------------ #
    # one sided: interested whether the random starts are less likely than the SO starts
    wilcoxon_rank_sums(background_transcripts_df, df_merged, 'best')

    # ------------------ PLOT SECOND BEST PROBABILITIES ------------------ #
    plot_SO_background(df_merged, background_transcripts_df,
                       'second_best', datatype, outdir)

    # ------------------ WILOCXON RANK SUM TEST FOR THE BEST STARTS ------------------ #
    wilcoxon_rank_sums(background_transcripts_df, df_merged, 'second_best')

    #################################################################################
    # ------------------ COMPARE WITH RIBOSEQ EMPIRICAL FINDINGS ------------------ #
    #################################################################################

    for empirical_Ribo_findings_file in glob.glob("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample/NMD_genome/*_unique_regions.csv"):
        Ribo_df_merged_first, Ribo_df_merged_second, sample = analyze_emp_background_Riboseq(
            empirical_Ribo_findings_file, df_merged)

        plot_emp_background_TransAI(
            Ribo_df_merged_first, Ribo_df_merged_second, outdir, datatype, sample)

    #################################################################################
    # ------------------ COMPARE WITH RIBOTISH VAL ORFS FINDINGS ------------------ #
    #################################################################################
    # create unique genomic region bedtool
    UR_BedTool = BedTool(
        f'/projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_{datatype}_16_12_24_chrom_sorted.bed')

    # get TranslationAI information required
    df_merged_RiboTISH = df_merged[[
        'TIS_dict', 'OrfTransID', 'TIS_pos_list', 'OrfStarts_sorted']]
    df_merged_RiboTISH = df_merged_RiboTISH.rename(
        columns={'OrfTransID': 'Tid_RT'})

    # Reformat to bed format, pybedtools
    for RiboTISH_file in glob.glob("/projects/splitorfs/work/Riboseq/Output/RiboTISH_NMD_custom/*.csv"):
        sample = os.path.basename(
            RiboTISH_file).rsplit('_', 3)[0]
        print(sample)
        RiboTISH_SO_df = pd.read_csv(RiboTISH_file, header=0)
        if len(RiboTISH_SO_df.index) > 0:
            RiboTISH_BedTool, RiboTISH_SO_df = create_RiboTISH_BedTool(
                RiboTISH_SO_df)

            URs_found_in_RiboTISH_df = obtain_correct_ORF_RiboTISH(
                UR_BedTool, RiboTISH_BedTool)

            _, URs_RiboTISH_TransAI_df_first, URs_RiboTISH_TransAI_df_second = get_ORF_start_probs(
                URs_found_in_RiboTISH_df, df_merged_RiboTISH, on='Tid_RT')

            plot_RiboTISH_TransAI(URs_RiboTISH_TransAI_df_first,
                                  URs_RiboTISH_TransAI_df_second, datatype, outdir, sample)


if __name__ == "__main__":
    main()
