# From a gffcompare tmap file, get the isoforms that have the same EJ chain
# as present in the reference "=" code, write them as txt

# usage: python get_equal_ejc_isoforms.py

import os
import os.path
import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get same EJ chain Isoforms of new assembly compared to ref by gffcompare.")

    parser.add_argument(
        "tmap_file",
        type=str,
        help="Path to tmap file of gffcompare"
    )

    return parser.parse_args()


def main(tmap_file):
    dir_name = os.path.dirname(tmap_file)
    sample_name = os.path.basename(tmap_file).split('.')[1]

    tmap_df = pd.read_csv(tmap_file, sep='\t')
    tmap_same_ejc_df = tmap_df[tmap_df['class_code'] == '=']
    tmap_same_ejc_df['qry_id'].to_csv(os.path.join(
        dir_name, f'{sample_name}_same_ejc_as_ref_isoforms.txt'), index=False)
    tmap_novel_df = tmap_df[tmap_df['class_code'] != '=']
    tmap_novel_df['qry_id'].to_csv(os.path.join(
        dir_name, f'{sample_name}_novel_isoforms.txt'), index=False)


if __name__ == "__main__":
    args = parse_arguments()
    main(args.tmap_file)
