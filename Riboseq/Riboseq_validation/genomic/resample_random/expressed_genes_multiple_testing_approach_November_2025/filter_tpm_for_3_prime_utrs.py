import os
import pandas as pd
import argparse
from rnanorm import FPKM, TPM
from sklearn import set_config
import seaborn as sbn
import matplotlib.pyplot as plt
import numpy as np
# set_config(transform_output='pandas')


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('three_prime_coverage_file',
                        help='Path to three prime coverage file')

    parser.add_argument('sample',
                        help='sample')

    parser.add_argument('tpm_threshold',
                        help='TPM Threshold')

    return parser.parse_args()


def main(three_prime_coverage_file, sample, tpm_threshold):
    outdir = os.path.dirname(three_prime_coverage_file)

    three_prime_coverage_file_df = pd.read_csv(
        three_prime_coverage_file, sep='\t', header=None,
        names=['chr', 'start', 'stop', '3_prime_name', 'score', 'strand', 'nr_overlap_reads', 'nr_bases_covered', 'length_3_prime_UTR', 'covered_fraction'])

    # TPM calculation in the 3' UTRs
    three_prime_coverage_file_df['RPK'] = three_prime_coverage_file_df['nr_overlap_reads'] / \
        three_prime_coverage_file_df['length_3_prime_UTR']

    scaling_factor = sum(three_prime_coverage_file_df['RPK'])/1000000

    three_prime_coverage_file_df['TPM'] = three_prime_coverage_file_df['RPK']/scaling_factor

    # plot TPM distribution
    ax = sbn.histplot(
        three_prime_coverage_file_df, x='TPM', bins=5000)
    ax.set_xlim(0, 400)
    ax.set_title(f'TPM values per 3 prime UTR for {sample}')
    plt.show()

    fig = ax.get_figure()  # get the figure that contains the axes
    fig.savefig(os.path.join(outdir,
                f'{sample}_TPM_0_400_3primes.png'), dpi=300, bbox_inches='tight')  # save as PNG

    print(f'Number of transcripts that have a TPM >= 100 for sample {sample}', sum(
        three_prime_coverage_file_df['TPM'] >= 100))

    # filter the regions for the genes with TPM above TPM Treshold
    three_primes_filtered = three_prime_coverage_file_df[three_prime_coverage_file_df['TPM']
                                                         <= tpm_threshold]

    three_primes_filtered.iloc[:, 0:6].to_csv(os.path.join(
        outdir, f'3_primes_tpm_filtered_{sample}.bed'), sep='\t', header=False, index=False)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    three_prime_coverage_file = args.three_prime_coverage_file
    sample = args.sample
    tpm_threshold = int(args.tpm_threshold)

    # three_prime_coverage_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/ERR3367797/3_primes_genomic_merged_numbered_ERR3367797_coverage.tsv'

    main(three_prime_coverage_file, sample, tpm_threshold)
