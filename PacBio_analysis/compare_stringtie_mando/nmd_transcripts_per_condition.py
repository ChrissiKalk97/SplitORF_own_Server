# ----- This script adds the assembly source to a GTF from TAMA merge ----- #
# ----- to the ninth column of the GTF with gene_id; transcript_id etc the ----- #
# ----- the original assembly will be added with source "assembly" ----- #

import os
import pandas as pd
import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("fifty_nt_csv",
                        help="Path to fifty_nt_csv")

    parser.add_argument("idxstats_dir",
                        help="Path to idxstats_dir")

    return parser.parse_args()


def main(fifty_nt_csv, idxstats_dir):

    fifty_nt_csv_df = pd.read_csv(fifty_nt_csv, index_col=0)
    nmd_transcript_list = fifty_nt_csv_df[fifty_nt_csv_df['50_nt'] == 1].index.to_list(
    )

    nmd_stats_per_sample_dict = {}
    for file in os.listdir(idxstats_dir):
        if file.endswith('idxstats.out'):
            idxpath = os.path.join(idxstats_dir, file)
            sample = file.split('_')[0] + '_' + file.split('_')[1]
            print(sample)
            idx_df = pd.read_csv(idxpath, sep='\t', header=None, names=[
                                 'tID', 'length', 'mapped', 'unmapped'])
            idx_df = idx_df[idx_df['mapped'] > 0]
            print(f'number of transcript models to which {sample} has a primary alignment', len(
                idx_df.index))
            idx_nmd_df = idx_df[idx_df['tID'].isin(nmd_transcript_list)]
            idx_nmd_df.to_csv(os.path.join(
                idxstats_dir, f'{sample}_NMD_idxstats.csv'))
            nmd_stats_per_sample_dict[sample] = [len(idx_df.index), len(
                idx_nmd_df), len(idx_nmd_df)/len(idx_df.index)]

    sample_nmd_summary_df = pd.DataFrame.from_dict(nmd_stats_per_sample_dict, orient='index', columns=[
                                                   'nr_supported_transcript_models', 'nr_supported_nmd_models', 'ratio_nmd_all_transcript_models'])
    cell_type = sample.split('_')[0]
    sample_nmd_summary_df.to_csv(os.path.join(
        idxstats_dir, f'{cell_type}_NMD_summary_idxstats.csv'))

    # save the list of NMD transcripts as .txt
    with open(os.path.join(idxstats_dir, f'{cell_type}_all_fifty_nt_transcripts_tama.txt'), 'w') as f:
        for trans in nmd_transcript_list:
            f.write(f"{trans}\n")


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    fifty_nt_csv = args.fifty_nt_csv
    idxstats_dir = args.idxstats_dir

    main(fifty_nt_csv, idxstats_dir)
