
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("cds_ur_intersection_bed",
                        help="Path to Bedtool intersection file of unqiue regions with CDSs")

    return parser.parse_args()


def main(cds_ur_intersection_bed):
    # read file with CDS cooridnates
    ur_cds_intersection_df = pd.read_csv(
        cds_ur_intersection_bed, sep='\t', header=None, names=['chr', 'start', 'stop', 'UR_name', 'score', 'strand', 'CDS_chr', 'CDS_start', 'CDS_stop', 'CDS_name', 'CDS_score', 'CDSstrand', 'num_bp'])

    print('Number of unique regions that do overlap with a CDS',
          len(ur_cds_intersection_df['UR_name'].unique()))


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    cds_ur_intersection_bed = args.cds_ur_intersection_bed

    # cds_ur_intersection_bed='/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-11.30.56_NMD_transcripts_correct_TSL_ref/Unique_DNA_Regions_genomic_CDS_intersection.bed'
    main(cds_ur_intersection_bed)
