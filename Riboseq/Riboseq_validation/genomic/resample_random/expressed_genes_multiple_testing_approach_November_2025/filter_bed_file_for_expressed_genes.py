import os
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("bed_file_to_filter",
                        help="Path to BED file to filter")

    parser.add_argument("htseq_counts",
                        help="Path to HTSeq-counts")

    return parser.parse_args()


def main(bed_file_to_filter, htseq_counts):
    outdir = os.path.dirname(htseq_counts)
    print(outdir)
    sample = os.path.basename(htseq_counts).rstrip('NMD_htseq_counts.tsv')
    bed_file_name = os.path.basename(bed_file_to_filter).rstrip('.bed')

    bed_file_to_filter_df = pd.read_csv(
        bed_file_to_filter, sep='\t', header=None)
    bed_file_to_filter_df['Gene_name'] = bed_file_to_filter_df.iloc[:, 3].apply(
        lambda x: x.split('|')[0])

    htseq_counts_df = pd.read_csv(
        htseq_counts, sep='\t', header=None, names=['Gene_name', 'count'])
    # filter out for the non-gene columns
    htseq_counts_df = htseq_counts_df[~htseq_counts_df['Gene_name'].apply(
        lambda x: x.startswith('__'))]
    htseq_counts_df_filtered = htseq_counts_df[htseq_counts_df['count'] > 2]

    bed_file_filtered = bed_file_to_filter_df[bed_file_to_filter_df['Gene_name'].isin(
        htseq_counts_df_filtered['Gene_name'])]

    bed_file_filtered.iloc[:, 0:6].to_csv(os.path.join(
        outdir, f'{bed_file_name}_{sample}.bed'), sep='\t', header=False, index=False)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    bed_file_to_filter = args.bed_file_to_filter
    htseq_counts = args.htseq_counts

    # bed_file_to_filter = '/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed'

    # htseq_counts = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/huvec_dnor_2/huvec_dnor_2_NMD_htseq_counts.tsv'

    main(bed_file_to_filter, htseq_counts)
