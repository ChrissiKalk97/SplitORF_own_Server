# ----- This script changes the order of transcript_id and gene_id in the 9th ----- #
# ----- field of a GTF that was generated using SQANTI3 to gene_id; transcript_id  ----- #

import os
import pandas as pd
import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("path_to_sqanti_gtf",
                        help="Path to SQANTI GTF to be changed")

    parser.add_argument("output_gtf",
                        help="Path to output GTF")

    return parser.parse_args()


def main(path_to_gtf, output_gtf):
    def add_exon_number(row):
        if row[2] == 'exon':
            row[8] = row[8] + f' exon_number "{int(row["exon_number"])}";'
        return row[8]

    gtf_df = pd.read_csv(path_to_gtf, header=None, sep='\t')
    gtf_df.iloc[:, 8] = gtf_df.iloc[:, 8].apply(
        lambda x: x.split(';')[1] + '; ' + x.split(';')[0] + ';')
    gtf_df["exon_number"] = (
        gtf_df[gtf_df[2] == "exon"]               # only exons
        .groupby(8)
        .cumcount() + 1                 # start numbering at 1
    ).astype(int)
    gtf_df.iloc[:, 8] = gtf_df.apply(lambda x: add_exon_number(x), axis=1)

    gtf_df.iloc[:, 0:9].to_csv(output_gtf, sep='\t', index=False, header=False,
                               quoting=csv.QUOTE_NONE, escapechar='\\')


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_gtf = args.path_to_sqanti_gtf
    output_gtf = args.output_gtf

    main(path_to_gtf, output_gtf)
