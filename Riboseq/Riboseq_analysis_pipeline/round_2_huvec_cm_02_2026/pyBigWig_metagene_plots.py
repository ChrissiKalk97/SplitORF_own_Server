# average the Riboseq coverage of a ceratain condition over replicates
# length normalize transcript coverage and plot a metagene plot across
# all transcripts
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


# path_to_bw_files = '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files'
# transcript_fai = '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai'
# mane_transcripts_cds_bed = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed'
# condition = 'DNOR'
# out_path = '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files/metagene_plots'
# color = '#1eb0e6'


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

    def calculate_std_sem(metagene_df, std_profile, n_reps):
        metagene_df['average_coverage_plus_std'] = metagene_df['average_coverage'] + \
            metagene_df['stddev_coverage']
        metagene_df['average_coverage_minus_std'] = metagene_df['average_coverage'] - \
            metagene_df['stddev_coverage']

        metagene_df['sem'] = std_profile / np.sqrt(n_reps)
        metagene_df['average_coverage_plus_sem'] = metagene_df['average_coverage'] + \
            metagene_df['sem']
        metagene_df['average_coverage_minus_sem'] = metagene_df['average_coverage'] - \
            metagene_df['sem']
        return metagene_df

    def plot_metagene(metagene_df, out_path, condition, color, region_type):
        plt.figure(figsize=(10, 5))

        # Plot the average line
        plt.plot(metagene_df[f'{region_type}_position'],
                 metagene_df['average_coverage'], label='Average Coverage', color=color)

        # Fill between +std and -std
        plt.fill_between(
            metagene_df[f'{region_type}_position'],
            metagene_df['average_coverage_minus_std'],
            metagene_df['average_coverage_plus_std'],
            color=color,
            alpha=0.3,
            label='±1 STD'
        )

        plt.xlim(min(metagene_df[f'{region_type}_position']),
                 max(metagene_df[f'{region_type}_position']))
        plt.xlabel(f'{region_type} position (%)')
        plt.ylabel(f'Avg Coverage (CPM)')
        plt.title(f'Metagene plot {condition}')
        plt.legend()
        plt.tight_layout()

        plt.savefig(os.path.join(out_path, f'{condition}_{region_type}_metagene_plot.svg'),
                    format='svg', dpi=300, bbox_inches='tight')
        plt.close()

    def calculate_and_plot_enrichment_df(transcript_coords, bw_files, out_path, condition, color, region_type):
        def normalize_transcript(signal, bins=100):
            """
            Normalize a transcript or region to a fixed number of bins.

            Parameters:
                signal : list or array
                    Coverage values along the transcript
                bins : int
                    Number of bins to scale to (default 100)

            Returns:
                np.ndarray of length `bins` with mean coverage per bin
            """
            signal = np.asarray(signal)

            if len(signal) >= bins:
                # If long enough, split into equal bins and average
                split = np.array_split(signal, bins)
                return np.array([chunk.mean() for chunk in split])
            else:
                # If shorter than bins, interpolate to scale up
                old_indices = np.linspace(0, 1, len(signal))
                new_indices = np.linspace(0, 1, bins)
                return np.interp(new_indices, old_indices, signal)

        # aggregate length normalized CPM signal for all transcripts per replicate
        metagene_cov_list = []
        for bw_object in bw_files:
            observed_values = []
            for transcript in transcript_coords.index:
                start = transcript_coords.loc[transcript,
                                              f'{region_type}_start']
                length = transcript_coords.loc[transcript,
                                               f'{region_type}_stop']
                bw_values = bw_object.values(
                    transcript, start, length)

                if len(bw_values) > 0 and sum(bw_values) > 5:
                    norm_transcript_cov = normalize_transcript(bw_values)
                    observed_values.append(norm_transcript_cov)

                else:
                    pass
                    # print(transcript, bw_object)

            # average all transcripts of one replicate
            try:
                if len(observed_values) > 1:
                    observed_values_np = np.array(observed_values)
                    observed_values_averaged_np = np.mean(
                        observed_values_np, axis=0)
                    metagene_cov_list.append(observed_values_averaged_np)

            except:
                pass

        # average replicates
        if len(metagene_cov_list) > 0:
            print(len(metagene_cov_list))
            all_reps = np.vstack(metagene_cov_list)
            mean_profile = np.mean(all_reps, axis=0)
            std_profile = np.std(all_reps, axis=0, ddof=1)
            # metaplot = np.mean(metagene_cov_list, axis=0)

            # assert len(bw_values_smoothed) == length
            metagene_df = pd.DataFrame({f'{region_type}_position': range(0, 100),
                                        'average_coverage': mean_profile,
                                        'stddev_coverage': std_profile
                                        })

            metagene_df = calculate_std_sem(
                metagene_df, std_profile, all_reps.shape[0])

            metagene_df.to_csv(os.path.join(
                path_to_bw_files, 'coordinates_per_transcript_csvs', f'{condition}_coords_for_plotting.csv'))
            plot_metagene(metagene_df, out_path, condition, color, region_type)

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
