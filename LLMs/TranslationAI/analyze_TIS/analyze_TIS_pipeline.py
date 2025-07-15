# ------------------ IMPORTS ------------------ #
import os
import argparse


from data_loader import load_TranslationAI, load_so_results, load_dna_ur_df
from helper_functions_analysis import nr_trans_and_mean_probs, get_TransAI_so_info, merge_df, \
    preprocess_k4neo_data, get_trans_ai_so_preds_df, \
    explode_so_df, subset_UR_for_expressed_genes, subset_validated_sos_df, val_so_by_position, all_URs_by_position, \
    val_perc_first_middle_last_orfs_csv, get_so_position_in_transcript, \
    compare_so_set_probabilities_by_position, identify_overlapping_unique_regions
from analysis_steps import perform_so_background_analysis, validated_so_per_sample_analysis, \
    RiboTISH_analysis
from plotting import plot_start_prob_by_orf_position, \
    plot_val_so_sets


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process Ribo-seq and TIS data for analysis."
    )

    # Required positional arguments
    parser.add_argument("TIS_results", help="Path to TIS results file")
    parser.add_argument("so_results", help="Path to SO results file")
    parser.add_argument("Ribo_coverage_path",
                        help="Path to Ribo-seq coverage directory")
    parser.add_argument("RiboTISH_path", help="Path to RiboTISH output")

    # Optional arguments
    parser.add_argument('--genes_to_keep_file', default="",
                        help='Path to the txt file with gene IDs to consider')
    parser.add_argument("--k4neo_path", default="",
                        help="Optional path to K4Neo results")
    parser.add_argument("--Ensembl_canonical_path", default="",
                        help="Optional path to Ensembl canonical file")

    return parser.parse_args()


def main(TIS_results, so_results, Ribo_coverage_path, RiboTISH_path, region_type, Ensembl_canonical_path='', k4neo_path='', genes_to_keep_file=''):
    UR_path = f'/projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_{region_type}_16_12_24_chrom_sorted.bed'

    # ------------------ DATA IMPORT ------------------ #
    TIS_results_df = load_TranslationAI(TIS_results)
    predicted_so_orfs, so_transcripts = load_so_results(so_results)

    # ------------------ COMPARE PREDICTIONS SO AND TRANSLATIONAI ------------------ #
    TransAI_so_preds = TIS_results_df[TIS_results_df['OrfTransID'].isin(
        so_transcripts)].copy()

    nr_trans_and_mean_probs(TransAI_so_preds, TIS_results_df, so_transcripts)

    TransAI_so_preds = get_TransAI_so_info(TransAI_so_preds)

    df_merged = merge_df(TransAI_so_preds, predicted_so_orfs)

    # ----------- COMPARE SO TO BACKGROUND TRANSLATIONAI PREDICTIONS --------- #
    perform_so_background_analysis(
        TIS_results_df, so_transcripts, df_merged, region_type, os.path.join(outdir, 'plots'))

    all_predicted_so_orfs, predicted_so_orfs, total_nr_so = explode_so_df(
        predicted_so_orfs)

    #################################################################################
    # ------------------ COMPARE WITH RIBOSEQ EMPIRICAL FINDINGS ------------------ #
    #################################################################################
    all_predicted_so_orfs, Ensembl_canonical_df = validated_so_per_sample_analysis(
        df_merged, Ribo_coverage_path, os.path.join(outdir, 'plots', 'Riboseq_single_samples'), region_type, all_predicted_so_orfs, Ensembl_canonical_path=Ensembl_canonical_path)

    validated_so_df, nr_validated_so, nr_validated_transcripts = subset_validated_sos_df(
        all_predicted_so_orfs)

    # ------------------ LOAD DNA UNIQUE REGIONS ------------------ #
    dna_ur_df, nr_orfs_with_UR, nr_transcripts_with_UR = load_dna_ur_df(
        UR_path)

    category_names = ['# validated SO',
                      '# SO with UR not validated', '# SO without UR']
    # ------------------ LOAD DNA UNIQUE REGIONS ------------------ #
    if len(genes_to_keep_file) > 0:
        nr_orfs_with_UR, dna_ur_df = subset_UR_for_expressed_genes(
            genes_to_keep_file, dna_ur_df, validated_so_df, df_merged)
        # define category_names
        # define nr_orfs_with_UR as nr_orfs_with_UR minus the ones that are not expressed
        # define
        category_names = [
            '# validated SO', '# SO with UR not validated but gene expressed', '# SO without UR or gene not expressed']
        region_type = f'{region_type}_gene_filtered'

    # ------------------ PLOT VALIDATED TRANSCRIPT PROPORTIONS ------------------ #
    plot_val_so_sets(nr_orfs_with_UR, nr_validated_so,
                     total_nr_so, os.path.join(outdir, 'plots'), region_type, category_names=category_names)

    # ------------------ PLOT VALIDATED SO POSITIONS ------------------ #
    nr_val_first_orfs, nr_val_middle_orfs, nr_val_last_orfs = val_so_by_position(
        validated_so_df, nr_validated_so, os.path.join(outdir, 'plots'), region_type)

    # ------------------ DNA URs among first, middle and last ORFs ------------------ #
    nr_first_orfs_ur, nr_middle_orfs_ur, nr_last_orfs_ur = all_URs_by_position(
        dna_ur_df, all_predicted_so_orfs, os.path.join(outdir, 'plots'), region_type)

    val_perc_first_middle_last_orfs_csv(nr_val_first_orfs,
                                        nr_val_middle_orfs,
                                        nr_val_last_orfs,
                                        nr_first_orfs_ur,
                                        nr_middle_orfs_ur,
                                        nr_last_orfs_ur,
                                        os.path.join(outdir, 'plots'),
                                        'validated_percentages_by_position.csv')

    # ------------------ TransSI probs of first, middle, last orfs ------------------ #
    val_so_trans_ai_df = get_trans_ai_so_preds_df(
        TransAI_so_preds, validated_so_df)
    all_predicted_so_orfs = get_so_position_in_transcript(
        all_predicted_so_orfs)
    all_so_trans_ai_df = get_trans_ai_so_preds_df(
        TransAI_so_preds, all_predicted_so_orfs)

    if len(genes_to_keep_file) > 0:
        compare_so_set_probabilities_by_position(
            all_so_trans_ai_df, dna_ur_df, val_so_trans_ai_df, region_type, os.path.join(outdir, 'plots'), gene_filtering=True)
    else:
        compare_so_set_probabilities_by_position(
            all_so_trans_ai_df, dna_ur_df, val_so_trans_ai_df, region_type, os.path.join(outdir, 'plots'))

    val_dna_overlapping_ur_df = identify_overlapping_unique_regions(
        validated_so_df, dna_ur_df, os.path.join(outdir, 'files'))

    plot_start_prob_by_orf_position(val_so_trans_ai_df,
                                    'Probability of Riboseq validated ORFs by ORF position',
                                    region_type,
                                    os.path.join(outdir, 'plots'))

    #################################################################################
    # ------------------ COMPARE WITH RIBOTISH VAL ORFS FINDINGS ------------------ #
    #################################################################################
    # create unique genomic region bedtool
    UR_BedTool = RiboTISH_analysis(
        df_merged, region_type, RiboTISH_path, os.path.join(outdir, 'plots'), UR_path)

    #################################################################################
    # ------------------ COMPARE WITH K4NEO TRANSCRIPTS          ------------------ #
    #################################################################################

    if k4neo_path:
        # k4neo was from an old SO pipeline run, where many more (false) unique regions were
        # found as the longest pc isoforms were used and not all pc isoforms with good tsl
        k4neo_TransAI_preds, k4neo_validated_transcripts = preprocess_k4neo_data(
            k4neo_path, UR_BedTool, df_merged)

        perform_so_background_analysis(
            TransAI_so_preds, k4neo_validated_transcripts, k4neo_TransAI_preds, region_type, os.path.join(outdir, 'plots'), so='k4neo')


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    TIS_results = args.TIS_results
    so_results = args.so_results
    Ribo_coverage_path = args.Ribo_coverage_path
    RiboTISH_path = args.RiboTISH_path
    genes_to_keep_file = args.genes_to_keep_file
    k4neo_path = args.k4neo_path
    Ensembl_canonical_path = args.Ensembl_canonical_path

    # TIS_results = '/projects/splitorfs/work/LLMs/TranslationAI/Output/NMD_trnascripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt'
    # so_results = '/projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt'
    # Ribo_coverage_path = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/NMD_genome'
    # genes_to_keep_file = '/projects/splitorfs/work/LLMs/TranslationAI/Output/analyze_TIS/genes_to_keep/genes_to_keep.txt'
    # RiboTISH_path = '/projects/splitorfs/work/Riboseq/Output/RiboTISH_NMD_custom'
    # Ensembl_canonical_path = '/projects/splitorfs/work/LLMs/TranslationAI/Input/NMD_transcripts_cDNA_coordinates.txt'
    # k4neo_path = '/projects/splitorfs/work/LLMs/TranslationAI/Input/k4neo_val_transcripts/NMD_transcripts_found_with_k4neo.txt'

    resultdir = os.path.dirname(TIS_results)
    region_type = os.path.basename(so_results).split('_')[1].split('.')[0]

    os.makedirs(f'{resultdir}/analyze_TIS', exist_ok=True)
    os.makedirs(f'{resultdir}/analyze_TIS/{region_type}', exist_ok=True)
    os.makedirs(f'{resultdir}/analyze_TIS/{region_type}/plots', exist_ok=True)
    os.makedirs(
        f'{resultdir}/analyze_TIS/{region_type}/plots/Riboseq_single_samples', exist_ok=True)
    os.makedirs(f'{resultdir}/analyze_TIS/{region_type}/files', exist_ok=True)
    outdir = f'{resultdir}/analyze_TIS/{region_type}'
    os.makedirs(f'{outdir}/plots/RiboTISH', exist_ok=True)

    main(TIS_results, so_results, Ribo_coverage_path,
         RiboTISH_path, region_type, Ensembl_canonical_path, k4neo_path, genes_to_keep_file)
