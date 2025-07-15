import os
import glob

import pandas as pd

from data_loader import load_Ensembl_canonical

from helper_functions_analysis import background_probs, wilcoxon_rank_sums, \
    analyze_emp_background_riboseq, preprocess_RiboTISH_files, \
    create_RiboTISH_BedTool, obtain_correct_orf_RiboTISH,  get_orf_start_probs, \
    check_start_probs_Ensembl_canonical

from plotting import plot_so_background, plot_validated_so_probs_TransAI, plot_RiboTISH_TransAI, \
    plot_RiboTISH_inframecount_vs_probs


def perform_so_background_analysis(TIS_results_df, so_transcripts, df_merged, datatype, outdir, so='so'):
    '''
    compares the best and second best start TransAI probabilities between 2 sets
    of transcripts
    '''
    background_transcripts_df = background_probs(
        TIS_results_df, so_transcripts)

    # ------------------ PLOT BEST PROBABILITIES ------------------ #
    plot_so_background(df_merged, background_transcripts_df,
                       'best', datatype, outdir, so)

    # ------------------ WILOCXON RANK SUM TEST FOR THE BEST STARTS ------------------ #
    # one sided: interested whether the random starts are less likely than the so starts
    wilcoxon_rank_sums(background_transcripts_df, df_merged, 'best')

    # ------------------ PLOT SECOND BEST PROBABILITIES ------------------ #
    plot_so_background(df_merged, background_transcripts_df,
                       'second_best', datatype, outdir, so)

    # ------------------ WILOCXON RANK SUM TEST FOR THE BEST STARTS ------------------ #
    wilcoxon_rank_sums(background_transcripts_df, df_merged, 'second_best')


def validated_so_per_sample_analysis(df_merged, Ribo_coverage_path, outdir, datatype, all_predicted_so_orfs, Ensembl_canonical_path=''):
    '''
    for coverage based validated SOs perform background analysis
    '''
    for empirical_Ribo_findings_file in glob.glob(f"{Ribo_coverage_path}/*_unique_regions.csv"):
        Ribo_df_merged, Ribo_df_merged_first, Ribo_df_merged_second, sample = analyze_emp_background_riboseq(
            empirical_Ribo_findings_file, df_merged)

        plot_validated_so_probs_TransAI(
            Ribo_df_merged_first, Ribo_df_merged_second, outdir, datatype, sample)

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

        Ensembl_canonical_df = None
        if Ensembl_canonical_path:
            Ensembl_canonical_df = load_Ensembl_canonical(
                Ensembl_canonical_path)
            check_start_probs_Ensembl_canonical(
                Ribo_df_merged, Ensembl_canonical_df, verbose=True)
    return all_predicted_so_orfs, Ensembl_canonical_df


def RiboTISH_analysis(df_merged, datatype, RiboTISH_path, outdir, UR_path):
    UR_BedTool, df_merged_RiboTISH = preprocess_RiboTISH_files(
        df_merged, datatype, UR_path)

    # Reformat to bed format, pybedtools
    for RiboTISH_file in glob.glob(f"{RiboTISH_path}/*.csv"):
        sample = os.path.basename(
            RiboTISH_file).rsplit('_', 3)[0]
        RiboTISH_so_df = pd.read_csv(RiboTISH_file, header=0)
        if len(RiboTISH_so_df.index) > 0:
            RiboTISH_BedTool, RiboTISH_so_df = create_RiboTISH_BedTool(
                RiboTISH_so_df)

            URs_found_in_RiboTISH_df = obtain_correct_orf_RiboTISH(
                UR_BedTool, RiboTISH_BedTool)

            RiboTISH_TransAI_df, URs_RiboTISH_TransAI_df_first, URs_RiboTISH_TransAI_df_second = get_orf_start_probs(
                URs_found_in_RiboTISH_df, df_merged_RiboTISH, on='Tid_RT')

            plot_RiboTISH_TransAI(URs_RiboTISH_TransAI_df_first,
                                  URs_RiboTISH_TransAI_df_second, datatype, outdir, sample)

            plot_RiboTISH_inframecount_vs_probs(
                RiboTISH_TransAI_df, datatype, outdir, sample)

    return UR_BedTool
