# calculate coefficient of variation for each transcript
# average across samples, report median CV and plot a histogram
import os
import argparse
import re

from collections import defaultdict


import pyBigWig
import numpy as np
import pandas as pd

import seaborn as sbn
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('--path_to_bw_files', help='Path to BigWig files')
    parser.add_argument(
        '--mane_gtf', help='Path to BED file of MANE CDS coords')
    parser.add_argument('--out_path', help='Path where to output the plots')
    parser.add_argument('--condition', help='condition to be plotted')
    parser.add_argument('--color', help='color of plot')
    parser.add_argument('--region_type', help='cds or transcript')

    return parser.parse_args()


# path_to_bw_files = '/projects/splitorfs/work/UPF1_deletion/Output/alignment_genome/STAR/deduplicated/coverage_bw_files'
# mane_gtf = '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf'
# condition = '0h'
# out_path = '/projects/splitorfs/work/UPF1_deletion/Output/alignment_genome/STAR/deduplicated/coverage_bw_files/cv_plots'
# color = '#1eb0e6'


def main(path_to_bw_files,  mane_gtf,
         out_path, condition, color='blue', region_type='transcript'):

    def get_bw_file_paths(path_to_bw_files, condition):
        bw_list = []
        for file in os.listdir(path_to_bw_files):
            if file.endswith('.bw') and condition in file:
                bw_list.append(
                    os.path.join(path_to_bw_files, file))
        return bw_list

    def read_in_bw(importin_over_input_bw_list):
        bw_files = []
        for bw_path in importin_over_input_bw_list:
            bw_object = pyBigWig.open(os.path.join(
                bw_path))
            bw_files.append(bw_object)
        return bw_files

    def get_coords_transcripts_of_interest(mane_gtf):
        # Example: loading GTF to get CDS coordinates per transcript
        gtf_df = pd.read_csv(mane_gtf, sep='\t', comment='#', header=None,
                             names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

        # filter for CDS
        cds_df = gtf_df[gtf_df['feature'] == 'CDS']
        cds_df['transcript_id'] = cds_df['attribute'].str.extract(
            'transcript_id "([^"]+)"')

        # parse attributes to get transcript_id
        cds_df['transcript_id'] = cds_df['attribute'].str.extract(
            'transcript_id "([^"]+)"')
        return cds_df

    def get_transcript_intervals(cds_df):
        # I do not care about the orientation, as the CV is only about neighboring bases
        transcript_cds_coords = defaultdict(list)
        for _, row in cds_df.iterrows():
            transcript_cds_coords[row['transcript_id']].append(
                (row['chrom'], row['start'], row['end'], row['strand'])
            )
        return transcript_cds_coords

    def plot_histogram(transcript_coords, out_path, condition, color, region_type):
        fig, ax = plt.subplots(figsize=(10, 5))
        # plt.figure(figsize=(10, 5))
        sbn.histplot(transcript_coords[transcript_coords > 0]['cv'], ax=ax)

        ax.set_title(f'coefficient of variation in transcripts in {condition}')
        ax.set_xlabel('coefficient of variation')

        plt.xlim(0, 8)

        plt.savefig(os.path.join(out_path, f'{condition}_{region_type}_cv_hist.svg'),
                    format='svg', dpi=300, bbox_inches='tight')
        plt.close()

    def calculate_and_plot_enrichment_df(transcript_cds_coords, bw_files, out_path, condition, color, region_type):
        def calculate_cv(signal):
            """calculates coefficient of variation for supplied signal as np.array
            Args:
                signal (np.arrray): CPM normalized signal per base
            Returns:
                cv: coefficient of variation for the respective transcript
            """
            signal = np.asarray(signal)
            if signal.mean() > 0.0:
                cv = np.std(signal) / np.mean(signal)
            return cv

        transcripts_cv_df = pd.DataFrame({'tid': transcript_cds_coords.keys()})
        transcripts_cv_df['cv'] = 0
        transcripts_cv_df = transcripts_cv_df.set_index('tid')
        for tx, intervals in transcript_cds_coords.items():
            observed_values = []
            for bw_object in bw_files:
                signal_values = []
                interval_counter = 0
                for chrom, start, end, strand in intervals:
                    # need to remove 'chr' from chromosome name: mapping Ensembl transcriptome
                    chrom = chrom[3:]
                    # BigWig is 0-based, end-exclusive
                    try:
                        if interval_counter == 0:
                            # start-1 if your GTF is 1-based
                            vals = bw_object.values(chrom, start - 21, end)
                        elif interval_counter == len(intervals) - 1:
                            vals = bw_object.values(chrom, start - 1, end + 20)
                        else:
                            vals = bw_object.values(chrom, start - 1, end)
                        # convert None to 0 (pyBigWig may return None for missing regions)
                        vals = [v if v is not None else 0.0 for v in vals]
                        signal_values.extend(vals)
                        interval_counter += 1
                    except RuntimeError:
                        try:
                            vals = bw_object.values(chrom, start-1, end)
                            interval_counter += 1
                        except:
                            continue
                if len(signal_values) > 0 and sum(signal_values) > 5:
                    cv_transcript = calculate_cv(signal_values)
                    observed_values.append(cv_transcript)

                else:
                    observed_values.append(0)
                    # print(transcript, bw_object)

                # average each transcript across replicates
                try:
                    if len(observed_values) > 1:
                        observed_values_np = np.array(observed_values)
                        observed_values_averaged_np = np.mean(
                            observed_values_np, axis=0)
                        transcripts_cv_df.loc[tx,
                                              'cv'] = observed_values_averaged_np
                except:
                    pass

        os.makedirs(os.path.join(
            path_to_bw_files, 'trans_cv_dfs'), exist_ok=True)
        os.makedirs(os.path.join(
            path_to_bw_files, 'cv_plots'), exist_ok=True)
        transcripts_cv_df.to_csv(os.path.join(
            path_to_bw_files, 'trans_cv_dfs', f'{condition}_trans_cv_df.csv'))
        print(
            f'Median coefficient of variation among transcripts for {condition}:', transcripts_cv_df['cv'].median())
        print(f'Nr of transcripts for {condition} with CPM > 5:', sum(
            transcripts_cv_df['cv'] > 0))
        plot_histogram(transcripts_cv_df, out_path,
                       condition, color, region_type)

    bw_list = get_bw_file_paths(
        path_to_bw_files, condition)

    bw_files = read_in_bw(bw_list)

    cds_df = get_coords_transcripts_of_interest(
        mane_gtf)

    transcript_cds_coords = get_transcript_intervals(cds_df)

    calculate_and_plot_enrichment_df(
        transcript_cds_coords, bw_files, out_path, condition, color, region_type)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_bw_files = args.path_to_bw_files
    mane_gtf = args.mane_gtf
    out_path = args.out_path
    condition = args.condition
    color = args.color
    region_type = args.region_type

    main(path_to_bw_files, mane_gtf,
         out_path, condition, color, region_type)
