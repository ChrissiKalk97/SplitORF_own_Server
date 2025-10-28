import os
import argparse

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("path_to_coordinate_csv",
                        help="Path to coordinate CSV files")
    parser.add_argument("--enrichment_length",
                        help="Nr of bp enrichment needs to persist")
    parser.add_argument("--enrichment_threshold",
                        help="Threshold enrichment needs to surpass")
    return parser.parse_args()


# path_to_coordinate_csv = "/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/bowtie1/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig/coordinates_per_transcript_csvs"
# # # coordinates are in codons
# enrichment_length = 50//3
# enrichment_threshold = 2


def main(path_to_coordinate_csv, enrichment_length, enrichment_threshold):
    onset_dict = {}
    for coordinate_csv in os.listdir(path_to_coordinate_csv):
        if coordinate_csv.endswith('.csv'):
            coordinate_df = pd.read_csv(os.path.join(
                path_to_coordinate_csv, coordinate_csv))
            coordinate_above_threshold_df = coordinate_df[coordinate_df['average_ratio']
                                                          > enrichment_threshold]
            if len(coordinate_above_threshold_df.index) > 0:
                nr_consecutive = 0
                consecutive = False
                consecutive_start = coordinate_above_threshold_df.iloc[0, 1]
                i = 0
                while consecutive == False and i < len(coordinate_above_threshold_df.index) - 1:
                    i += 1
                    if coordinate_above_threshold_df.iloc[i, 1] == coordinate_above_threshold_df.iloc[i-1, 1] + 1:
                        nr_consecutive += 1
                    else:
                        nr_consecutive = 0
                        consecutive_start = coordinate_above_threshold_df.iloc[i + 1, 1]
                    if nr_consecutive == enrichment_length:
                        consecutive = True
                        onset_dict[coordinate_csv] = [consecutive_start]
            if coordinate_csv not in onset_dict:
                onset_dict[coordinate_csv] = [None]
    onset_df = pd.DataFrame.from_dict(onset_dict).transpose(copy=True)
    onset_df['name'] = onset_df.index
    onset_df = onset_df.sort_values(by=['name'], axis=0)
    onset_df.to_csv(
        os.path.join(path_to_coordinate_csv, 'onset_per_csv.csv'), index=False)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_coordinate_csv = args.path_to_coordinate_csv
    enrichment_length = int(args.enrichment_length)
    enrichment_threshold = int(args.enrichment_threshold)

    main(path_to_coordinate_csv, enrichment_length, enrichment_threshold)
