# ------------------ IMPORTS ------------------ #
import os
import os.path
import argparse


from data_loader import load_so_results, load_dna_ur_df
from helper_functions_analysis import explode_so_df, \
    subset_UR_for_expressed_genes, subset_validated_sos_df, \
    val_so_by_position, all_URs_by_position, \
    val_perc_first_middle_last_orfs_csv, count_orfs_by_position, \
    identify_overlapping_unique_regions, validated_so_per_sample_analysis, \
    add_sample_info_ur_df
from plotting import plot_val_so_sets, plot_three_category_pie


def parse_args():
    parser = argparse.ArgumentParser(
        description="Form validated set of URs from Ribo-seq data taking into account positioning of the URs."
    )

    # Required positional arguments
    parser.add_argument("--so_results", help="Path to SO results file")
    parser.add_argument("--ribo_coverage_path",
                        help="Path to Ribo-seq coverage directory")
    parser.add_argument("--ur_path",
                        help="Path to Unique region genomic BED file")
    parser.add_argument("--result_dir",
                        help="Directory for results")
    parser.add_argument("--region_type",
                        help="Directory for results")

    return parser.parse_args()


def main(so_results, ribo_coverage_path, region_type, ur_path, outdir):
    # ------------------ DATA IMPORT ------------------ #

    predicted_so_orfs, so_transcripts = load_so_results(so_results)

    all_predicted_so_orfs, predicted_so_orfs, total_nr_so = explode_so_df(
        predicted_so_orfs)

    #################################################################################
    # ------------------ COMPARE WITH RIBOSEQ EMPIRICAL FINDINGS ------------------ #
    #################################################################################
    all_predicted_so_orfs = validated_so_per_sample_analysis(
        ribo_coverage_path, all_predicted_so_orfs)

    validated_so_df, nr_validated_so, nr_validated_transcripts = subset_validated_sos_df(
        all_predicted_so_orfs, outdir, region_type)

    # ------------------ LOAD DNA UNIQUE REGIONS ------------------ #
    dna_ur_df, nr_orfs_with_UR, nr_transcripts_with_UR = load_dna_ur_df(
        ur_path)

    # ------------------ LOAD DNA UNIQUE REGIONS ------------------ #
    nr_orfs_with_UR, dna_ur_df = subset_UR_for_expressed_genes(
        dna_ur_df, validated_so_df, ribo_coverage_path, outdir, region_type)

    # define category_names
    # define nr_orfs_with_UR as nr_orfs_with_UR minus the ones that are not expressed
    # define
    category_names = [
        '# validated SO', '# SO with UR not validated but gene expressed', '# SO without UR or gene not expressed']
    region_type = f'{region_type}'

    # ------------------ PLOT VALIDATED TRANSCRIPT PROPORTIONS ------------------ #
    plot_val_so_sets(nr_orfs_with_UR, nr_validated_so,
                     total_nr_so, os.path.join(outdir, 'plots'), region_type, category_names=category_names)

    print(f'Validation percentage for {region_type}:',
          nr_validated_so/nr_orfs_with_UR)

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
                                        outdir,
                                        'validated_percentages_by_position.csv')

    val_dna_overlapping_ur_df = identify_overlapping_unique_regions(
        validated_so_df, dna_ur_df, outdir)

    # need to combine both dfs: val_dna_overlapping_ur_df for UR overlap info and
    # validated_so_df for sample information

    total_nr_val_distinct_orfs = len(val_dna_overlapping_ur_df.index)

    nr_first_orfs, nr_middle_orfs, nr_last_orfs = count_orfs_by_position(
        val_dna_overlapping_ur_df)
    plot_three_category_pie(nr_first_orfs,
                            nr_middle_orfs,
                            nr_last_orfs,
                            total_nr_val_distinct_orfs,
                            ['# first ORFs', '# middle ORFs', '# last ORFs'],
                            'Pie chart of ORF positions of SOs with DNA UR',
                            os.path.join(outdir, 'plots'),
                            'positions_of_distinct_UR_validated_pie_chart',
                            region_type,
                            ['#75C1C5', '#FFC500', '#CC79A7']
                            )

    val_dna_overlapping_ur_df = add_sample_info_ur_df(
        validated_so_df, val_dna_overlapping_ur_df, outdir)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    so_results = args.so_results
    ribo_coverage_path = args.ribo_coverage_path
    ur_path = args.ur_path
    result_dir = args.result_dir
    region_type = args.region_type

    # so_results = '/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-13.07.17_NMD_for_paper/UniqueProteinORFPairs.txt'
    # ribo_coverage_path = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome'
    # ur_path = '/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-13.07.17_NMD_for_paper/Unique_DNA_Regions_genomic_final.bed'
    # result_dir = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis'

    os.makedirs(f'{result_dir}/SO_valdiation', exist_ok=True)
    os.makedirs(f'{result_dir}/SO_valdiation/{region_type}', exist_ok=True)
    os.makedirs(
        f'{result_dir}/SO_valdiation/{region_type}/plots', exist_ok=True)
    outdir = f'{result_dir}/SO_valdiation/{region_type}'

    main(so_results, ribo_coverage_path, region_type, ur_path, outdir)
