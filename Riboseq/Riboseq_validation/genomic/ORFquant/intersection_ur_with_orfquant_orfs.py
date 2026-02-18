import pandas as pd
import os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=".")

    parser.add_argument(
        "intersection_bed",
        type=str,
        help="Path to intersection bed of URs with ORF GTF regions from ORFquant"
    )
    return parser.parse_args()


def main(intersection_bed):
    intersection_df = pd.read_csv(intersection_bed, sep='\t', header=None, names=[
        'UR chr',
        'UR start',
        'UR end',
        'UR name',
        'Ur score',
        'UR strand',
        'ORF chr',
        'ORF start',
        'ORF end',
        'ORf gene name',
        'ORf score',
        'ORF strand',
        'ORfQuant',
        'exon',
        'another score',
        'transcript_information',
        'nr_bp_overlap'])

    intersection_df = intersection_df[[
        'UR chr',
        'UR start',
        'UR end',
        'UR name',
        'UR strand',
        'ORF chr',
        'ORF start',
        'ORF end',
        'ORf gene name',
        'ORF strand',
        'transcript_information',
        'nr_bp_overlap'
    ]].copy()

    intersection_df = intersection_df[intersection_df['UR strand']
                                      == intersection_df['ORF strand']]
    intersection_df['fraction_overlap'] = intersection_df.loc[:, 'nr_bp_overlap'] / \
        (intersection_df.loc[:, 'UR end'] - intersection_df.loc[:, 'UR start'])

    intersection_df = intersection_df[intersection_df['fraction_overlap'] == 1.0].copy(
    )
    # 1900 URs validated...das glaubt doch kein Mensch....
    sample = os.path.basename(intersection_bed)
    sample_grouped = sample[:-3] + '_completely_overlapping_grouped.csv'
    sample = sample[:-3] + '_completely_overlapping.csv'
    outdir = os.path.dirname(intersection_bed)
    intersection_df.to_csv(os.path.join(outdir, sample))

    intersection_df_grouped = intersection_df.groupby('UR name').agg({'UR chr': 'first',
                                                                      'UR start': 'first',
                                                                      'UR end': 'first',

                                                                      'UR strand': 'first',
                                                                      'ORF chr': 'first',
                                                                      'ORF start': list,
                                                                      'ORF end': list,
                                                                      'ORf gene name': list,
                                                                      'ORF strand': 'first',
                                                                      'transcript_information': list,
                                                                      'nr_bp_overlap': list}).reset_index()
    # 876 still sounds like very much

    intersection_df_grouped.to_csv(os.path.join(outdir, sample_grouped))


if __name__ == "__main__":
    args = parse_arguments()
    main(args.intersection_bed)
