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
    parser.add_argument("path_to_csv_files",
                        help="Path to CSV files per transcript and comparison")
    parser.add_argument("transcript_to_plot",
                        help="the name of the transcript that should be plotted")
    parser.add_argument("out_path", help="Path where to output the plots")
    parser.add_argument("importin", help="importin to be plotted")
    parser.add_argument(
        "--puro", help="Set Puro to 'Puro' if Puro samples are to be plotted")
    parser.add_argument(
        "--background1", help="")
    parser.add_argument(
        "--background2", help="")
    parser.add_argument(
        "--color1", help="set color for the first csv file")
    parser.add_argument(
        "--color2", help="set color for the second csv file")

    return parser.parse_args()


# path_to_csv_files = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig/coordinates_per_transcript_csvs'
# transcript_to_plot = 'ENST00000397885.3'
# out_path = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript/plot_input_and_mock'
# importin = 'A2'
# background1 = 'In'
# background2 = 'M'
# puro = ''
# color1 = '#1eb0e6'
# color2 = 'orange'


def plot_enrichment(transcript_df, out_path, transcript, importin, background1, background2, color1, color2):
    plt.figure(figsize=(10, 5))

    # Plot the average line
    plt.plot(transcript_df['CDS_position'],
             transcript_df['average_ratio'], label=f'Average {importin} IP/{background1}', color=color1)

    # Fill between +std and -std
    plt.fill_between(
        transcript_df['CDS_position'],
        transcript_df['min_value'],
        transcript_df['max_value'],
        color=color1,
        alpha=0.3,
        label=f'Min Max value IP/{background1}'
    )

    # Plot the average line
    plt.plot(transcript_df['CDS_position'],
             transcript_df[f'average_ratio_{background2}'], label=f'Average {importin} IP/{background2}', color=color2)

    # Fill between +std and -std
    plt.fill_between(
        transcript_df['CDS_position'],
        transcript_df[f'min_value_{background2}'],
        transcript_df[f'max_value_{background2}'],
        color=color2,
        alpha=0.3,
        label=f'Min Max value IP/{background2}'
    )

    plt.axhline(y=2.0, color='black', linestyle='--',
                linewidth=1, label='FC Threshold = 2')
    plt.ylim(0, 20)
    plt.xlim(min(transcript_df['CDS_position']),
             max(transcript_df['CDS_position']))
    plt.xlabel('Transcript position (Codons)')
    plt.ylabel(f'Enrichment ({importin} IP/{background1} or {background2})')
    plt.title(f'SeRP enrichment signal across {transcript} CDS')
    plt.legend()
    plt.tight_layout()
    if puro == 'Puro':
        pass
    else:
        plt.savefig(os.path.join(out_path, f'{transcript}_enrichment_{importin}_over_{background1}_vs_{background2}_CDS_only_from_whole_trans.svg'),
                    format='svg', dpi=300, bbox_inches='tight')
    plt.close()


def main(path_to_csv_files, transcript_to_plot,
         out_path, importin,  puro='', background1='In', background2='M',
         color1='blue', color2='orange'):

    for file in os.listdir(path_to_csv_files):
        if file.endswith('.csv'):
            if file.startswith(transcript_to_plot) and importin in file and background1 in file:
                enrichment1_df = pd.read_csv(os.path.join(
                    path_to_csv_files, file), header=0, index_col=0)
            elif file.startswith(transcript_to_plot) and importin in file and background2 in file:
                enrichment2_df = pd.read_csv(os.path.join(
                    path_to_csv_files, file), header=0, index_col=0)

    enrichment1_df = enrichment1_df[[
        'CDS_position', 'average_ratio', 'min_value', 'max_value']].copy()

    enrichment1_df[f'average_ratio_{background2}'] = enrichment2_df['average_ratio']
    enrichment1_df[f'min_value_{background2}'] = enrichment2_df['min_value']
    enrichment1_df[f'max_value_{background2}'] = enrichment2_df['max_value']

    plot_enrichment(enrichment1_df, out_path, transcript_to_plot,
                    importin, background1, background2, color1, color2)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_csv_files = args.path_to_csv_files
    transcript_to_plot = args.transcript_to_plot
    out_path = args.out_path
    importin = args.importin
    puro = args.puro
    background1 = args.background1
    background2 = args.background2
    color1 = args.color1
    color2 = args.color2

    main(path_to_csv_files, transcript_to_plot,
         out_path, importin, puro, background1, background2, color1, color2)
