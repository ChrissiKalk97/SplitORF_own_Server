# This script takes the validated unique regions and checks whether their genes are
# RBPs, it write a TXT file of validated RBP genes as well as the validated_so_df
# filtered for RBPs and with the RBP information
# ------------------ IMPORTS ------------------ #
import os
import os.path
import argparse
import glob
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description='Summarize the trans.'
    )

    # Required positional arguments

    parser.add_argument('--rbp_file',
                        help='Path to RNA-binding protein CSV from RBPDB')
    parser.add_argument('--validation_file',
                        help='of the validated_so_df.csv which is checked for RBP validations')
    parser.add_argument('--region_type',
                        help='NMD or RI')

    return parser.parse_args()


def main(rbp_file, validation_file, region_type):
    outdir = os.path.dirname(rbp_file)
    os.makedirs(os.path.join(outdir, region_type), exist_ok=True)

    rbp_df = pd.read_csv(rbp_file, sep=';')
    rbp_df = rbp_df.iloc[:, 1:].copy()
    validation_df = pd.read_csv(validation_file)
    validation_df = validation_df.iloc[:, 1:].copy()

    rbp_df_filtered = rbp_df[rbp_df['Annotation ID'].isin(
        validation_df['geneID'])].copy()
    print(f'Number of different RBPs validated {region_type}', len(rbp_df_filtered.index))

    rbp_df_filtered = rbp_df_filtered.rename(
        columns={'Annotation ID': 'geneID'})

    rbp_val_df = validation_df.merge(rbp_df_filtered, on='geneID', how='right')
    rbp_val_df.to_csv(os.path.join(outdir, region_type, 'rbp_val_df.csv'))

    rbp_val_df['geneID'].unique().tofile(os.path.join(
        outdir, region_type, 'validated_rbp_genes.txt'), sep='\n')


if __name__ == '__main__':
    args = parse_args()

    rbp_file = args.rbp_file
    validation_file = args.validation_file
    region_type = args.region_type

    # rbp_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/proteins.php'
    # validation_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/validated_so_df.csv'
    # region_type = 'NMD'

    main(rbp_file, validation_file, region_type)
