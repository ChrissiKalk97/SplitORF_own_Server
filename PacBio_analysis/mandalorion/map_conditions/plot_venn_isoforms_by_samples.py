from upsetplot import UpSet
from upsetplot import from_contents
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import os


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Summarize the number of mapped reads for each isoform per sample using idxstats.")

    parser.add_argument(
        "isoform_by_sample_csv_path",
        type=str,
        help="Path to the idxstats files"
    )
    return parser.parse_args()


def main(isoform_by_sample_csv_path):
    isoform_by_sample_df = pd.read_csv(isoform_by_sample_csv_path, index_col=0)
    isoform_by_sample_df = isoform_by_sample_df.iloc[:-1, :]
    upset_dict = {}
    for sample_col in isoform_by_sample_df.iloc[:, 2:].columns:
        sample_subset_df = isoform_by_sample_df.loc[:, [
            "isoform", sample_col]].copy()
        sample_subset_df = sample_subset_df[sample_subset_df[sample_col] > 2].copy(
        )
        upset_dict[sample_col] = sample_subset_df['isoform'].to_list()

    isoforms_by_sample_upset = from_contents(upset_dict)
    ax_dict = UpSet(isoforms_by_sample_upset, subset_size="count").plot()
    outname = os.path.basename(isoform_by_sample_csv_path)
    outname = outname.split('_')[0]
    outdir = os.path.dirname(isoform_by_sample_csv_path)
    plt.savefig(os.path.join(
        outdir, f'{outname}_isoform_by_sample_upset.svg'), format="svg")


if __name__ == "__main__":
    args = parse_arguments()
    main(args.isoform_by_sample_csv_path)
