# From a gffcompare tmap file, get the isoforms that have the same EJ chain
# as present in the reference "=" code, write them as txt

# usage: python get_equal_ejc_isoforms.py

import os
import os.path
import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get the number of novel NMD transcripts from mandalorion assembly.")

    parser.add_argument(
        "fifty_nt_csv",
        type=str,
        help="Path to CSV of fifty_nt pipe run with mando assembly"
    )

    parser.add_argument(
        "gffcompare_novel_transcripts",
        type=str,
        help="Path to TXT file containing novel (not '=') transcripts in mando ass determined by gffcompare"
    )

    parser.add_argument(
        "--assembly_type", default='',
        type=str,
        help="Indicate assembly comapred to, e.g. full or filtered"
    )

    return parser.parse_args()


def main(fifty_nt_csv, gffcompare_novel_transcripts, assembly_type):
    out_dir = os.path.dirname(fifty_nt_csv)
    sample_name = os.path.basename(fifty_nt_csv).split('.')[0]

    fifty_nt_df = pd.read_csv(fifty_nt_csv, index_col=0)
    nmd_positive_df = fifty_nt_df[fifty_nt_df['50_nt'] == 1.0]

    novel_transcripts_df = pd.read_csv(gffcompare_novel_transcripts)
    nmd_positive_df.index.isin(novel_transcripts_df['qry_id'])
    novel_nmd_isoforms_df = nmd_positive_df[nmd_positive_df.index.isin(
        novel_transcripts_df['qry_id'])]
    novel_nmd_isoforms_df.to_csv(os.path.join(
        out_dir, f'{sample_name}_{assembly_type}_novel_nmd_isoforms.txt'), columns=[], header=False)


if __name__ == "__main__":
    args = parse_arguments()
    main(args.fifty_nt_csv, args.gffcompare_novel_transcripts, args.assembly_type)
