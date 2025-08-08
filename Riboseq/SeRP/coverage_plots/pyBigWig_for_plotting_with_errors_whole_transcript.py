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
    parser.add_argument("mane_transcripts_cds_bed",
                        help="Path to CDS coordinates BED file of MANE transcripts")
    parser.add_argument("transcripts_to_plot_txt",
                        help="Path to TXT file with transcripts to plot")
    parser.add_argument("out_path", help="Path where to output the plots")
    parser.add_argument("importin", help="importin to be plotted")
    parser.add_argument("background", help="background to be plotted")
    parser.add_argument(
        "--puro", help="Set Puro to 'Puro' if Puro samples are to be plotted")
    parser.add_argument(
        "--color", help="")

    return parser.parse_args()


path_to_bw_files = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig'
mane_transcripts_cds_bed = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed'
transcripts_to_plot_txt = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt'
out_path = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript'
importin = 'A2'
background = 'In'
puro = ''
color = '#1eb0e6'


def main(path_to_bw_files, mane_transcripts_cds_bed, transcripts_to_plot_txt,
         out_path, importin='A2', background='Input', puro='', color='blue'):

    def get_bw_file_paths(path_to_bw_files, importin, background, puro):
        if puro == 'Puro':
            importin_over_input_bw_list = []
            for file in os.listdir(path_to_bw_files):
                if file.endswith('.bw'):
                    if file.startswith(importin) and background in file and 'Puro' in file:
                        importin_over_input_bw_list.append(
                            os.path.join(path_to_bw_files, file))
        else:
            importin_over_input_bw_list = []
            for file in os.listdir(path_to_bw_files):
                if file.endswith('.bw'):
                    if file.startswith(importin) and background in file and 'CHX' in file:
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

    def get_coords_transcripts_of_interest(transcripts_to_plot_txt, mane_transcripts_cds_bed):
        transcripts_of_interest_df = pd.read_csv(
            transcripts_to_plot_txt, header=None, names=['tid'])
        MANE_CDS_coords = pd.read_csv(
            mane_transcripts_cds_bed, sep='\t', header=None, names=['tid', 'start', 'stop'])
        MANE_CDS_coords = MANE_CDS_coords[MANE_CDS_coords['tid'].isin(
            transcripts_of_interest_df['tid'])]
        MANE_CDS_coords = MANE_CDS_coords.set_index('tid')
        return MANE_CDS_coords

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

    def plot_enrichment(transcript_df, out_path, transcript, importin, background, color):
        plt.figure(figsize=(10, 5))

        # Plot the average line
        plt.plot(transcript_df['CDS_position'],
                 transcript_df['average_ratio'], label='Average Ratio', color=color)

        # Fill between +std and -std
        plt.fill_between(
            transcript_df['CDS_position'],
            transcript_df['average_ratio_minus_sem'],
            transcript_df['average_ratio_plus_sem'],
            color=color,
            alpha=0.3,
            label='±1 SEM'
        )

        plt.axhline(y=2.0, color='black', linestyle='--',
                    linewidth=1, label='FC Threshold = 2')
        plt.ylim(0, 20)
        plt.xlim(min(transcript_df['CDS_position']),
                 max(transcript_df['CDS_position']))
        plt.xlabel('Transcript position (Codons)')
        plt.ylabel(f'Enrichment (IP {importin}/{background})')
        plt.title(f'SeRP enrichment signal across {transcript} CDS')
        plt.legend()
        plt.tight_layout()
        if puro == 'Puro':
            plt.savefig(os.path.join(out_path, f'{transcript}_enrichment_{importin}_over_{background}_{puro}_CDS_only_from_whole_trans.svg'),
                        format='svg', dpi=300, bbox_inches='tight')
        else:
            plt.savefig(os.path.join(out_path, f'{transcript}_enrichment_{importin}_over_{background}_CDS_only_from_whole_trans.svg'),
                        format='svg', dpi=300, bbox_inches='tight')
        plt.close()

    def calculate_and_plot_enrichment_df(MANE_CDS_coords, bws_to_average, out_path, importin, background, color):
        for transcript in MANE_CDS_coords.index:
            start = MANE_CDS_coords.loc[transcript, 'start']
            stop = MANE_CDS_coords.loc[transcript, 'stop']
            observed_values = []
            for bw_object in bws_to_average:
                bw_values = bw_object.values(transcript, start, stop)
                # can make this as a hyperparameter for the respective window size...
                bw_values = [1.0] * 9 + bw_values + [1.0] * 9
                # calculate the average of each three bp
                bw_values_binned = [sum(bw_values[i:i+3]) /
                                    3 for i in range(0, len(bw_values), 3)]
                # write the average for 3bp each: bw_values_binned[i//3] is the same number for always 3
                # number: floor division

                bw_values_smoothed = [
                    sum(bw_values_binned[i-3:i+4])/7 for i in range(3, len(bw_values_binned)-3, 1)]

                observed_values.append(bw_values_smoothed)

            assert len({len(values) for values in observed_values}) == 1
            observed_values_np = np.array(observed_values)
            observed_values_averaged_np = np.mean(observed_values_np, axis=0)
            observed_std_np = np.std(observed_values_np, axis=0)
            observed_values_min_np = np.min(observed_values_np, axis=0)
            observed_values_max_np = np.max(observed_values_np, axis=0)

            transcript_df = pd.DataFrame({'CDS_position': range(0, len(bw_values_smoothed)),
                                          'average_ratio': observed_values_averaged_np,
                                          'stddev_ratio': observed_std_np,
                                          'min_value': observed_values_min_np,
                                          'max_value': observed_values_max_np})

            transcript_df = calculate_std_sem(transcript_df, observed_values)
            if puro == 'Puro':
                transcript_df.to_csv(os.path.join(
                    path_to_bw_files, 'coordinates_per_transcript_csvs', f'{transcript}_{importin}_{background}_{puro}_coords_for_plotting.csv'))
            else:
                transcript_df.to_csv(os.path.join(
                    path_to_bw_files, 'coordinates_per_transcript_csvs', f'{transcript}_{importin}_{background}_coords_for_plotting.csv'))
            plot_enrichment(transcript_df, out_path,
                            transcript, importin, background, color)

    importin_over_input_bw_list = get_bw_file_paths(
        path_to_bw_files, importin, background, puro)

    bws_to_average = read_in_bw(importin_over_input_bw_list)

    MANE_CDS_coords = get_coords_transcripts_of_interest(
        transcripts_to_plot_txt, mane_transcripts_cds_bed)

    calculate_and_plot_enrichment_df(
        MANE_CDS_coords, bws_to_average, out_path, importin, background, color)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_bw_files = args.path_to_bw_files
    mane_transcripts_cds_bed = args.mane_transcripts_cds_bed
    transcripts_to_plot_txt = args.transcripts_to_plot_txt
    out_path = args.out_path
    importin = args.importin
    background = args.background
    puro = args.puro
    color = args.color

    main(path_to_bw_files, mane_transcripts_cds_bed, transcripts_to_plot_txt,
         out_path, importin, background, puro, color)
