# This script calculates a csv file with the nr of mapped reads per sample for each isoform
# the required input is the path to the idxstat files of the long reads aligned to the transcriptome

import pandas as pd
import os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Summarize the number of mapped reads for each isoform per sample using idxstats.")

    parser.add_argument(
        "idx_path",
        type=str,
        help="Path to the idxstats files"
    )
    return parser.parse_args()


def main(idx_path):
    idx_files = [os.path.join(idx_path, f) for f in os.listdir(
        idx_path) if f.endswith('idxstats.out')]

    idx_summary_df = pd.read_csv(idx_files[0], sep='\t', header=None)
    idx_summary_df = idx_summary_df.rename(
        columns={0: 'isoform', 1: 'length', 2: f'mapped', 3: 'unmapped'})
    idx_summary_df = idx_summary_df.loc[:, ['isoform', 'length']].copy()

    for nr, idx in enumerate(idx_files):
        sample = os.path.basename(idx)
        sample = sample.removesuffix("_pbmm2_aligned_idxstats.out")
        idx_df = pd.read_csv(idx, sep='\t', header=None)
        nr = nr + 1
        idx_df = idx_df.rename(
            columns={0: 'isoform', 1: 'length', 2: f'sample_{nr}', 3: 'unmapped'})
        idx_dict = dict(zip(idx_df['isoform'], idx_df[f'sample_{nr}']))
        idx_summary_df[f'sample_{nr}'] = idx_summary_df['isoform'].map(
            idx_dict)

    sample = sample.split('_')[0]
    idx_summary_df = idx_summary_df.drop('length', axis=1)

    idx_summary_df = idx_summary_df.rename(
        columns={idx_summary_df.columns[0]: "superPBID"}).copy()
    idx_summary_df.to_csv(os.path.join(
        idx_path, f'{sample}_idx_fl_counts.txt'), sep='\t', index=False)


if __name__ == "__main__":
    args = parse_arguments()
    main(args.idx_path)
