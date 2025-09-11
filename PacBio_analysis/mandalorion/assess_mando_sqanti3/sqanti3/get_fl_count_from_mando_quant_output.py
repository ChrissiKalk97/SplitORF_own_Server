import os
import pandas as pd
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get fl count from Mando quant file.")

    parser.add_argument(
        "quant_path",
        type=str,
        help="Path to the input Mando GTF file"
    )

    return parser.parse_args()


# huvec_path = "/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/HUVEC"
# cm_path = "/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/CM"


def main(quant_path):
    quant_file = quant_path + '/Isoforms.filtered.clean.quant'
    quant_df = pd.read_csv(quant_file, sep='\s+', header=0)
    quant_df = quant_df.drop(columns=['Gene']).copy()
    quant_df = quant_df.rename(
        columns={samp_number: 'sample' + samp_number for samp_number in quant_df.columns})
    quant_df = quant_df.rename(columns={"sampleIsoform": "superPBID"})
    sample = os.path.basename(quant_path)
    quant_df.to_csv(os.path.join(
        quant_path, sample + '_fl_counts.tsv'), sep='\t', index=False, header=True)


if __name__ == "__main__":
    args = parse_arguments()
    main(args.quant_path)
