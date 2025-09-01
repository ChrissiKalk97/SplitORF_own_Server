import os
import pandas as pd
import seaborn as sbn
import argparse
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("path_to_refine_csv_summaries",
                        help="Path to BigWig files")

    return parser.parse_args()


def main(path_to_refine_csv_summaries):
    csv_summaries = [file for file in os.listdir(
        path_to_refine_csv_summaries) if '.csv' in file]
    for file in csv_summaries:
        sample_name = '_'.join(file.split('_')[0:2])

        csv_path = os.path.join(path_to_refine_csv_summaries, file)
        refine_df = pd.read_csv(csv_path)

        # Example: histogram of the "insertlen" column
        sbn.histplot(refine_df['insertlen'], bins=100,
                     kde=False, color="skyblue")

        plt.xlabel("Insert length (bp)")
        plt.ylabel("Count")
        plt.title(f"Distribution of FLNC read lengths in {sample_name}")
        plt.show()

        plt.savefig(os.path.join(path_to_refine_csv_summaries,
                    f"{sample_name}_insertlen_after_refine_hist.png"),
                    dpi=300, bbox_inches="tight")

        plt.close()


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_refine_csv_summaries = args.path_to_refine_csv_summaries

    main(path_to_refine_csv_summaries)
