# ------------------ IMPORTS ------------------ #
import sys
import os
import os.path
import argparse
import glob
import pandas as pd

# fmt: off
# caution: path[0] is reserved for script path (or '' in REPL)
# insert path to helper functions from SO val set analysis to reuse the functionality
sys.path.insert(1, '/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/expressed_genes_multiple_testing_approach_November_2025/SO_validated_set_analysis/')


from helper_functions_analysis import calculate_overlapping_region_percentage
# fmt: on


def parse_args():
    parser = argparse.ArgumentParser(
        description='Summarize the trans.'
    )

    # Required positional arguments

    parser.add_argument('--ribo_coverage_path',
                        help='Path to Ribo-seq coverage directory')
    parser.add_argument('--out_name',
                        help='Name of the output summary CSV file')
    parser.add_argument('--samples_of_interest_string',
                        help='String with unique parts of the samples to be \
                        included, separated by comma, e.g. \
                            "sample1,sample2"')

    return parser.parse_args()


def main(ribo_coverage_path, out_name, samples_of_interest_string):
    samples_of_interest = samples_of_interest_string.split(',')

    transcripts_with_two_urs_dfs = []
    for empirical_ribo_findings_file in glob.glob(f'{ribo_coverage_path}/*_two_regions_validated_on_transcript.csv'):
        if any(sample_name in empirical_ribo_findings_file for sample_name in samples_of_interest):
            sample = os.path.basename(
                empirical_ribo_findings_file).removesuffix('_two_regions_validated_on_transcript.csv')
            ribo_results_df = pd.read_csv(
                empirical_ribo_findings_file, header=0, index_col=0)
            ribo_results_df['sample'] = sample
            for transcript in ribo_results_df['name'].unique():
                ribo_results_df_transcript = ribo_results_df[ribo_results_df['name'] == transcript]
                start1 = ribo_results_df_transcript.iloc[0, 4]
                end1 = ribo_results_df_transcript.iloc[0, 5]
                start2 = ribo_results_df_transcript.iloc[1, 4]
                end2 = ribo_results_df_transcript.iloc[1, 5]
                ribo_results_df['ORF'] = ribo_results_df['new_name'].apply(
                    lambda x: x.split(':')[1])
                if len(ribo_results_df_transcript.index) == 2 and len(ribo_results_df['ORF'].unique()) > 1:
                    overlap = calculate_overlapping_region_percentage(
                        start1, end1, start2, end2)
                    if overlap == 0:
                        transcripts_with_two_urs_dfs.append(
                            ribo_results_df_transcript)
                elif len(ribo_results_df_transcript.index) > 2 and len(ribo_results_df_transcript.index) < 5 and len(ribo_results_df['ORF'].unique()) > 1:
                    start3 = ribo_results_df_transcript.iloc[2, 4]
                    end3 = ribo_results_df_transcript.iloc[2, 5]
                    overlap1 = calculate_overlapping_region_percentage(
                        start1, end1, start2, end2)
                    overlap2 = calculate_overlapping_region_percentage(
                        start1, end1, start3, end3)
                    overlap3 = calculate_overlapping_region_percentage(
                        start3, end3, start2, end2)
                    if overlap1 == 0 or overlap2 == 0 or overlap3 == 0:
                        transcripts_with_two_urs_dfs.append(
                            ribo_results_df_transcript)
                    if len(ribo_results_df_transcript.index) == 4:
                        start4 = ribo_results_df_transcript.iloc[3, 4]
                        end4 = ribo_results_df_transcript.iloc[3, 5]
                        overlap4 = calculate_overlapping_region_percentage(
                            start1, end1, start4, end4)
                        overlap5 = calculate_overlapping_region_percentage(
                            start4, end4, start2, end2)
                        overlap6 = calculate_overlapping_region_percentage(
                            start3, end3, start4, end4)
                        if overlap1 != 0 and overlap2 != 0 and overlap3 != 0 and (overlap4 == 0 or overlap5 == 0 or overlap6 == 0):
                            transcripts_with_two_urs_dfs.append(
                                ribo_results_df_transcript)
                elif len(ribo_results_df_transcript.index) > 4:
                    print('Transcript with > 4 URs validated',
                          ribo_results_df_transcript)

    summarized_df = pd.concat(transcripts_with_two_urs_dfs, ignore_index=True)
    agg_rule = {col: 'first' for col in summarized_df.columns}
    agg_rule['sample'] = lambda x: ','.join(x)
    summarized_df = summarized_df.groupby('new_name').agg(agg_rule)
    summarized_df.to_csv(out_name)


if __name__ == '__main__':
    args = parse_args()

    ribo_coverage_path = args.ribo_coverage_path
    out_name = args.out_name
    samples_of_interest_string = args.samples_of_interest_string

    # out_name = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/transcripts_with_2_regions_summarized.csv'
    # ribo_coverage_path = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome'
    # samples_of_interest_string = 'SRR10,SRR85,HCT'

    main(ribo_coverage_path, out_name, samples_of_interest_string)
