import os
import argparse

import pyBigWig
import numpy as np
import pandas as pd

import seaborn as sbn
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("path_to_bw_files", help="Path to BigWig files")
    parser.add_argument(
        "transcript_fai", help="Path to FASTA fai index for transcripts")
    parser.add_argument("transcripts_to_plot_txt",
                        help="Path to TXT file with transcripts to plot")
    parser.add_argument("out_path", help="Path where to output the plots")
    parser.add_argument("--numerator", help="numerator to be plotted")
    parser.add_argument("--background", help="background to be plotted")
    parser.add_argument(
        "--color", help="")
    parser.add_argument(
        "--window_size", help="window size must be mutliple of 3 and should be an odd number multiple of 3")

    return parser.parse_args()


# path_to_bw_files = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/enrichment_plots_whole_trans'
# transcript_fai = '/projects/serp/work/references/Chaetomium_thermophilum_longest_transcript.fasta.fai'
# transcripts_to_plot_txt = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0_not_IP_WT_vs_IN_0.5_andpadj_0.05.txt'
# out_path = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/enrichment_plots_whole_trans/S_E_over_S_In_DEG_plots'
# numerator = 'S_E'
# background = 'S_In'
# color = '#1eb0e6'
# window_size = 63


def main(path_to_bw_files, transcript_fai, transcripts_to_plot_txt,
         out_path, numerator='S_E', background='S_In', color='blue', window_size=21):

    def get_bw_file_paths(path_to_bw_files, numerator, background):
        importin_over_input_bw_list = []
        for file in os.listdir(path_to_bw_files):
            if file.endswith('.bw'):
                if file.startswith(numerator) and background in file:
                    importin_over_input_bw_list.append(
                        os.path.join(path_to_bw_files, file))
        return importin_over_input_bw_list

    def read_in_bw(importin_over_input_bw_list):
        bws_to_average = []
        for bw_path in importin_over_input_bw_list:
            bw_object = pyBigWig.open(os.path.join(
                bw_path))
            bws_to_average.append(bw_object)
        return bws_to_average

    def get_coords_transcripts_of_interest(transcripts_to_plot_txt, transcript_fai):
        transcripts_of_interest_df = pd.read_csv(
            transcripts_to_plot_txt, header=None, names=['tid'])
        transcript_coords = pd.read_csv(
            transcript_fai, sep='\t', header=None, names=['tid', 'length', 'offset', 'linebases', 'linewidth', 'qualoffset'])
        transcript_coords = transcript_coords[transcript_coords['tid'].isin(
            transcripts_of_interest_df['tid'])]
        transcript_coords['start'] = 0
        transcript_coords = transcript_coords.set_index('tid')
        transcript_coords = transcript_coords[['start', 'length']]
        return transcript_coords

    def calculate_std_sem(transcript_df, observed_values):
        transcript_df['average_ratio_plus_std'] = transcript_df['average_ratio'] + \
            transcript_df['stddev_ratio']
        transcript_df['average_ratio_minus_std'] = transcript_df['average_ratio'] - \
            transcript_df['stddev_ratio']

        n = len(observed_values)
        transcript_df['sem'] = transcript_df['stddev_ratio'] / (n ** 0.5)
        transcript_df['average_ratio_plus_sem'] = transcript_df['average_ratio'] + \
            transcript_df['sem']
        transcript_df['average_ratio_minus_sem'] = transcript_df['average_ratio'] - \
            transcript_df['sem']
        return transcript_df

    def plot_enrichment(transcript_df, out_path, transcript, numerator, background, color, window_size=21):
        plt.figure(figsize=(10, 5))

        # Plot the average line
        plt.plot(transcript_df['transcript_position'],
                 transcript_df['average_ratio'], label='Average Ratio', color=color)

        # Fill between +std and -std
        plt.fill_between(
            transcript_df['transcript_position'],
            transcript_df['average_ratio_minus_sem'],
            transcript_df['average_ratio_plus_sem'],
            color=color,
            alpha=0.3,
            label='±1 SEM'
        )

        plt.axhline(y=2.0, color='black', linestyle='--',
                    linewidth=1, label='FC Threshold = 2')
        plt.ylim(0, 20)
        plt.xlim(min(transcript_df['transcript_position']),
                 max(transcript_df['transcript_position']))
        plt.xlabel('Transcript position (Codons)')
        plt.ylabel(f'Enrichment (IP {numerator}/{background})')
        plt.title(f'SeRP enrichment signal across {transcript} transcript')
        plt.legend()
        plt.tight_layout()

        plt.savefig(os.path.join(out_path, f'{transcript}_enrichment_{numerator}_over_{background}_ws{window_size}.svg'),
                    format='svg', dpi=300, bbox_inches='tight')
        plt.close()

    def calculate_and_plot_enrichment_df(transcript_coords, bws_to_average, out_path, numerator, background, color, window_size=21):
        buffer_zone = int((window_size - 1)/2)

        # idea: instead of the MANE CDS coords do faidx to get transcript length
        for transcript in transcript_coords.index:
            start = transcript_coords.loc[transcript, 'start']
            length = transcript_coords.loc[transcript, 'length']
            observed_values = []

            for bw_object in bws_to_average:
                bw_values = bw_object.values(
                    transcript, start, length)

                if len(bw_values) > 0:

                    bw_values = [1.0] * buffer_zone + \
                        bw_values + [1.0] * buffer_zone

                    bw_values_smoothed = [
                        sum(bw_values[i-buffer_zone:i+buffer_zone])/window_size for i in range(buffer_zone, len(bw_values)-buffer_zone, 1)]

                    observed_values.append(bw_values_smoothed)

                else:
                    print(transcript, bw_object)

                # assert len({len(values) for values in observed_values}) == 1
                try:
                    observed_values_np = np.array(observed_values)
                    observed_values_averaged_np = np.mean(
                        observed_values_np, axis=0)
                    observed_std_np = np.std(observed_values_np, axis=0)
                    observed_values_min_np = np.min(observed_values_np, axis=0)
                    observed_values_max_np = np.max(observed_values_np, axis=0)
                except:
                    print(transcript)
                    print(start)
                    print(length)
                    print(observed_values)
                    print(numerator, background)

                if len(observed_values[0]) > 0:

                    transcript_df = pd.DataFrame({'transcript_position': range(0, len(bw_values_smoothed)),
                                                  'average_ratio': observed_values_averaged_np,
                                                  'stddev_ratio': observed_std_np,
                                                  'min_value': observed_values_min_np,
                                                  'max_value': observed_values_max_np})

                    transcript_df = calculate_std_sem(
                        transcript_df, observed_values)

                    transcript_df.to_csv(os.path.join(
                        path_to_bw_files, 'coordinates_per_transcript_csvs', f'{transcript}_{numerator}_{background}_coords_for_plotting.csv'))
                    plot_enrichment(transcript_df, out_path,
                                    transcript, numerator, background, color, window_size)

                else:
                    print('No observed values for:', transcript)

    print('Lets go!')

    numerator_over_input_bw_list = get_bw_file_paths(
        path_to_bw_files, numerator, background)

    bws_to_average = read_in_bw(numerator_over_input_bw_list)

    transcript_coords = get_coords_transcripts_of_interest(
        transcripts_to_plot_txt, transcript_fai)

    calculate_and_plot_enrichment_df(
        transcript_coords, bws_to_average, out_path, numerator, background, color, window_size)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_bw_files = args.path_to_bw_files
    transcript_fai = args.transcript_fai
    transcripts_to_plot_txt = args.transcripts_to_plot_txt
    out_path = args.out_path
    numerator = args.numerator
    background = args.background
    color = args.color
    window_size = int(args.window_size)

    main(path_to_bw_files, transcript_fai, transcripts_to_plot_txt,
         out_path, numerator, background, color, window_size)
