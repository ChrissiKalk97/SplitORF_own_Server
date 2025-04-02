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

from data_loader import load_TranslationAI, load_SO_results, load_Ensembl_canonical
from analysis import nr_trans_and_mean_probs, get_TransAI_SO_info, merge_df, \
    background_probs, wilcoxon_rank_sums, analyze_emp_background_Riboseq, \
    create_RiboTISH_BedTool, obtain_correct_ORF_RiboTISH, get_ORF_start_probs, \
    check_start_probs_Ensembl_canonical
from plotting import plot_SO_background, plot_emp_background_TransAI, plot_RiboTISH_TransAI,\
                        plot_RiboTISH_inframecount_vs_probs






def main(TIS_results, SO_results, Ribo_coverage_path, RiboTISH_path, Ensembl_canonical_path='', k4neo_path=''):
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
    print(Ribo_coverage_path)
    for empirical_Ribo_findings_file in glob.glob(f"{Ribo_coverage_path}/*_unique_regions.csv"):
        print(empirical_Ribo_findings_file)
        Ribo_df_merged , Ribo_df_merged_first, Ribo_df_merged_second, sample = analyze_emp_background_Riboseq(
            empirical_Ribo_findings_file, df_merged)

        plot_emp_background_TransAI(
            Ribo_df_merged_first, Ribo_df_merged_second, outdir, datatype, sample)
        
        if Ensembl_canonical_path:
            Ensembl_canonical_df = load_Ensembl_canonical(Ensembl_canonical_path)
            check_start_probs_Ensembl_canonical(Ribo_df_merged, Ensembl_canonical_df, verbose = True)

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
    for RiboTISH_file in glob.glob(f"{RiboTISH_path}/*.csv"):
        sample = os.path.basename(
            RiboTISH_file).rsplit('_', 3)[0]
        print(sample)
        RiboTISH_SO_df = pd.read_csv(RiboTISH_file, header=0)
        if len(RiboTISH_SO_df.index) > 0:
            RiboTISH_BedTool, RiboTISH_SO_df = create_RiboTISH_BedTool(
                RiboTISH_SO_df)

            URs_found_in_RiboTISH_df = obtain_correct_ORF_RiboTISH(
                UR_BedTool, RiboTISH_BedTool)

            RiboTISH_TransAI_df, URs_RiboTISH_TransAI_df_first, URs_RiboTISH_TransAI_df_second = get_ORF_start_probs(
                URs_found_in_RiboTISH_df, df_merged_RiboTISH, on='Tid_RT')

            plot_RiboTISH_TransAI(URs_RiboTISH_TransAI_df_first,
                                  URs_RiboTISH_TransAI_df_second, datatype, outdir, sample)
            
            plot_RiboTISH_inframecount_vs_probs(RiboTISH_TransAI_df, datatype, outdir, sample)


    #################################################################################
    # ------------------ COMPARE WITH K4NEO TRANSCRIPTS          ------------------ #
    #################################################################################
    
    # k4neo was from an old SO pipeline run, where many mroe (false) unique regions were
    # found as the longest pc isoforms were used and not all pc isoforms with good tsl
    k4neo_validated_transcripts = pd.read_csv(k4neo_path, index_col = None, header = None, names = ['OrfTransID'])
    k4neo_validated_transcripts = k4neo_validated_transcripts['OrfTransID'].to_list()

    # get URs
    UR_names = []
    for interval in UR_BedTool:
        UR_names.append(interval.name)

    UR_names = [re.split(r'[:|]', name)[1] for name in UR_names]

    # filter k4neo transcripts for transcripts that still have a UR predicted with the 
    # current pipeline
    k4neo_validated_transcripts = [trans for trans in k4neo_validated_transcripts if trans in UR_names]


    k4neo_TransAI_preds = df_merged[df_merged['OrfTransID'].isin(k4neo_validated_transcripts)]
    Non_k4neo_TransAI_preds = df_merged[~df_merged['OrfTransID'].isin(k4neo_validated_transcripts)]

    k4neo_TransAI_preds['best_start'] = k4neo_TransAI_preds['best_start_SO']
    k4neo_TransAI_preds['second_best_start'] = k4neo_TransAI_preds['second_best_start_SO']
    Non_k4neo_TransAI_preds['best_start'] = Non_k4neo_TransAI_preds['best_start_SO']
    Non_k4neo_TransAI_preds['second_best_start'] = Non_k4neo_TransAI_preds['second_best_start_SO']

    # ------------------ PLOT BEST PROBABILITIES ------------------ #
    plot_SO_background(k4neo_TransAI_preds, Non_k4neo_TransAI_preds,
                       'best', f'k4neo_{datatype}', outdir, SO = 'k4neo')

    # ------------------ WILOCXON RANK SUM TEST FOR THE BEST STARTS ------------------ #
    # one sided: interested whether the random starts are less likely than the SO starts
    wilcoxon_rank_sums(Non_k4neo_TransAI_preds, k4neo_TransAI_preds, 'best')

    # ------------------ PLOT SECOND BEST PROBABILITIES ------------------ #
    plot_SO_background(k4neo_TransAI_preds, Non_k4neo_TransAI_preds,
                       'second_best', f'k4neo_{datatype}', outdir, SO = 'k4neo')
    # ------------------ WILOCXON RANK SUM TEST FOR THE BEST STARTS ------------------ #
    wilcoxon_rank_sums(Non_k4neo_TransAI_preds, k4neo_TransAI_preds, 'second_best')


if __name__ == "__main__":
    print('script is running...')
    # ------------------ CONSTANTS ------------------ #
    TIS_results = sys.argv[1]
    SO_results = sys.argv[2]
    Ribo_coverage_path = sys.argv[3]
    RiboTISH_path = sys.argv[4]

    k4neo_path = ''
    if len(sys.argv) > 5: 
        k4neo_path = sys.argv[5]

    Ensembl_canonical_path = ''
    if len(sys.argv) > 6: 
        Ensembl_canonical_path = sys.argv[6]




    # TIS_results = '/projects/splitorfs/work/LLMs/TranslationAI/Output/NMD_trnascripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt'
    # SO_results = '/projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt'
    # Ribo_coverage_path = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample/NMD_genome'
    # RiboTISH_path = '/projects/splitorfs/work/Riboseq/Output/RiboTISH_NMD_custom'
    # Ensembl_canonical_path = '/projects/splitorfs/work/LLMs/TranslationAI/Input/NMD_transcripts_cDNA_coordinates.txt'
    # k4neo_path = '/projects/splitorfs/work/LLMs/TranslationAI/Input/k4neo_val_transcripts/NMD_transcripts_found_with_k4neo.txt'


    resultdir = os.path.dirname(TIS_results)
    datatype = os.path.basename(SO_results).split('_')[1].split('.')[0]

    os.makedirs(f'{resultdir}/plots', exist_ok=True)
    os.makedirs(f'{resultdir}/plots/{datatype}', exist_ok=True)
    outdir = f'{resultdir}/plots/{datatype}'
    os.makedirs(f'{outdir}/RiboTISH', exist_ok=True)
    print(TIS_results, SO_results, Ribo_coverage_path, RiboTISH_path, Ensembl_canonical_path, k4neo_path)
    main(TIS_results, SO_results, Ribo_coverage_path, RiboTISH_path, Ensembl_canonical_path, k4neo_path)
