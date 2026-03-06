# calculate coefficient of variation for each transcript
# average across samples, report median CV and plot a histogram
import os
import argparse

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
        '--transcript_fai', help='Path to FASTA fai index for transcripts')
    parser.add_argument(
        '--mane_transcripts_cds_bed', help='Path to BED file of MANE CDS coords')
    parser.add_argument('--out_path', help='Path where to output the plots')
    parser.add_argument('--condition', help='condition to be plotted')
    parser.add_argument('--color', help='color of plot')
    parser.add_argument('--region_type', help='cds or transcript')

    return parser.parse_args()


path_to_bw_files = '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files'
transcript_fai = '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai'
mane_transcripts_cds_bed = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed'
condition = 'DNOR'
out_path = '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files/metagene_plots'
color = '#1eb0e6'


def main(path_to_bw_files, transcript_fai, mane_transcripts_cds_bed,
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

    def get_coords_transcripts_of_interest(mane_transcripts_cds_bed, transcript_fai):
        mane_cds_coords = pd.read_csv(
            mane_transcripts_cds_bed, sep='\t', header=None, names=['tid', 'cds_start', 'cds_stop'])
        transcript_coords = pd.read_csv(
            transcript_fai, sep='\t', header=None, names=['tid', 'length', 'offset', 'linebases', 'linewidth', 'qualoffset'])
        mane_cds_coords['transcript_start'] = 0
        trans_len_dict = dict(
            zip(transcript_coords['tid'], transcript_coords['length']))
        mane_cds_coords['transcript_stop'] = mane_cds_coords['tid'].map(
            trans_len_dict)
        mane_cds_coords = mane_cds_coords.set_index('tid')
        return mane_cds_coords

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

    def calculate_and_plot_enrichment_df(transcript_coords, bw_files, out_path, condition, color, region_type):
        def calculate_cv(signal):
            """calculates coefficient of variation for supplied signal as np.array
            Args:
                signal (np.arrray): CPM normalized signal per base
            Returns:
                cv: coefficient of variation for the respective transcript
            """
            signal = np.asarray(signal)
            if signal .mean() > 0:
                cv = np.std(signal) / np.mean(signal)
            return cv

        # aggregate length normalized CPM signal for all transcripts per replicate
        transcript_coords['cv'] = 0
        for transcript in transcript_coords.index:
            observed_values = []
            for bw_object in bw_files:

                start = transcript_coords.loc[transcript,
                                              f'{region_type}_start']
                length = transcript_coords.loc[transcript,
                                               f'{region_type}_stop']
                if region_type == 'cds':
                    try:
                        bw_values = bw_object.values(
                            transcript, start - 20, length + 20)
                    except:
                        bw_values = bw_object.values(
                            transcript, start, length)

                if len(bw_values) > 0 and sum(bw_values) > 5:
                    cv_transcript = calculate_cv(bw_values)
                    observed_values.append(cv_transcript)

                else:
                    observed_values.append(0)
                    # print(transcript, bw_object)

            # average all transcripts of one replicate
            try:
                if len(observed_values) > 1:
                    observed_values_np = np.array(observed_values)
                    observed_values_averaged_np = np.mean(
                        observed_values_np, axis=0)
                    transcript_coords.loc[transcript,
                                          'cv'] = observed_values_averaged_np

            except:
                pass
        os.makedirs(os.path.join(
            path_to_bw_files, 'trans_cv_dfs'), exist_ok=True)
        os.makedirs(os.path.join(
            path_to_bw_files, 'cv_plots'), exist_ok=True)
        transcript_coords.to_csv(os.path.join(
            path_to_bw_files, 'trans_cv_dfs', f'{condition}_trans_cv_df.csv'))
        print(
            f'Median coefficient of variation among transcripts for {condition}:', transcript_coords['cv'].median())
        print(f'Nr of transcripts for {condition} with CPM > 5:', sum(
            transcript_coords['cv'] > 0))
        plot_histogram(transcript_coords, out_path,
                       condition, color, region_type)

    bw_list = get_bw_file_paths(
        path_to_bw_files, condition)

    bw_files = read_in_bw(bw_list)

    transcript_coords = get_coords_transcripts_of_interest(
        mane_transcripts_cds_bed, transcript_fai)

    calculate_and_plot_enrichment_df(
        transcript_coords, bw_files, out_path, condition, color, region_type)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_bw_files = args.path_to_bw_files
    transcript_fai = args.transcript_fai
    mane_transcripts_cds_bed = args.mane_transcripts_cds_bed
    out_path = args.out_path
    condition = args.condition
    color = args.color
    region_type = args.region_type

    main(path_to_bw_files, transcript_fai, mane_transcripts_cds_bed,
         out_path, condition, color, region_type)
