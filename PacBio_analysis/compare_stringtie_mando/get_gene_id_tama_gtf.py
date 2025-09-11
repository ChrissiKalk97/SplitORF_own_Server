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
    parser.add_argument("path_to_tama_gtf",
                        help="Path to TAMA GTF to be changed")

    parser.add_argument("path_to_classif",
                        help="Path to classification file SQANTI")

    parser.add_argument("output_gtf",
                        help="Path to output GTF")

    return parser.parse_args()


def main(path_to_tama_gtf, path_to_classif, output_gtf):
    def change_gene_id(row):
        gene_info_string = row[8]
        gene_info_list = gene_info_string.split('"')
        gene_info_list[1] = row['gene_id']
        gene_info_string = '"'.join(gene_info_list)
        return gene_info_string

    tama_gtf_df = pd.read_csv(path_to_tama_gtf, header=None, sep='\t')
    classif_df = pd.read_csv(path_to_classif, header=0, sep='\t')
    tama_gtf_df['gene_id'] = tama_gtf_df.iloc[:, 8].apply(
        lambda x: x.split('"')[1])
    classif_df['gene_id'] = classif_df['isoform'].apply(
        lambda x: x.split('.')[0])
    gene_id_dict = dict(
        zip(classif_df['gene_id'], classif_df['associated_gene']))

    tama_gtf_df['gene_id'] = tama_gtf_df['gene_id'].map(gene_id_dict)
    tama_gtf_df.iloc[:, 8] = tama_gtf_df.apply(
        lambda row: change_gene_id(row), axis=1)
    tama_gtf_df.iloc[:, 0:9].to_csv(output_gtf, sep='\t', index=False, header=False,
                                    quoting=csv.QUOTE_NONE, escapechar='\\')


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_tama_gtf = args.path_to_tama_gtf
    output_gtf = args.output_gtf
    path_to_classif = args.path_to_classif

    main(path_to_tama_gtf, path_to_classif, output_gtf)
