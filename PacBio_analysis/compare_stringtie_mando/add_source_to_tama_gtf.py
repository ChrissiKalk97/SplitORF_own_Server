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
    parser.add_argument("tama_gtf",
                        help="Path to TAMA GTF to be changed")

    parser.add_argument("tama_trans_info_file",
                        help="Path to TAMA trans information file")

    return parser.parse_args()


def main(tama_gtf, tama_trans_info_file):

    gtf_df = pd.read_csv(tama_gtf, header=None, sep='\t')
    tama_trans_df = pd.read_csv(tama_trans_info_file, header=0, sep='\t')

    outdir = os.path.dirname(tama_gtf)
    outname = os.path.basename(tama_gtf)[:-4]

    # get transcript_id column for gtf_df
    gtf_df['transcript_id'] = gtf_df.iloc[:, 8].apply(
        lambda x: x.split(";")[1] if len(x.split(";")) > 1 else None)
    gtf_df['transcript_id'] = gtf_df['transcript_id'].apply(
        lambda x: x.split('"')[1] if len(x.split('"')) > 1 else None)
    # set transcript_id column as index for tama_trans_df to use it for accessing columns
    tama_trans_df = tama_trans_df.set_index('transcript_id')

    # add 'sources' to ninth gtf field
    gtf_df['sources'] = gtf_df.apply(lambda row: tama_trans_df.loc[row['transcript_id'],
                                     'sources'] if row['transcript_id'] != None else None, axis=1)
    gtf_df.iloc[:, 8] = gtf_df.apply(lambda row: row[8] + ' sources "' +
                                     row['sources'] + '";' if row['sources'] != None else row[8], axis=1)

    # add 'all_sources_trans' to 9th GTF field
    gtf_df['all_source_trans'] = gtf_df.apply(
        lambda row: tama_trans_df.loc[row['transcript_id'], 'all_source_trans'] if row['transcript_id'] != None else None, axis=1)
    gtf_df.iloc[:, 8] = gtf_df.apply(lambda row: row[8] + ' all_source_trans "' +
                                     row['all_source_trans'] + '";' if row['all_source_trans'] != None else row[8], axis=1)

    gtf_df.iloc[:, 0:9].to_csv(os.path.join(outdir, f'{outname}.gtf'), sep='\t', index=False, header=False,
                               quoting=csv.QUOTE_NONE, escapechar='\\')


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    tama_gtf = args.tama_gtf
    tama_trans_info_file = args.tama_trans_info_file

    main(tama_gtf, tama_trans_info_file)
