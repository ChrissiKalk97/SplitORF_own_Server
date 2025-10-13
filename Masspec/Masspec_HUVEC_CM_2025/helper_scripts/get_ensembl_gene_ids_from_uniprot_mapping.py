import os
import pandas as pd
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Uniprot_Ensembl_mapping_file_tsv")

    parser.add_argument(
        "--uniprot_ensembl_tsv",
        required=True,
        help="."
    )
    return parser.parse_args()


# uniprot_ensembl_tsv = '/Users/christina/Documents/own_data/Masspec/SO_with_reference_masspec_files_16_09_25/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/Uniprot_to_ensembl_all_genes_HUVEC.tsv'


def main(uniprot_ensembl_tsv):
    out_dir = os.path.dirname(uniprot_ensembl_tsv)
    uniprot_ensembl_df = pd.read_csv(uniprot_ensembl_tsv, sep='\t')
    uniprot_ensembl_df['To'] = uniprot_ensembl_df['To'].apply(
        lambda x: x.split('.')[0])
    uniprot_ensembl_unique_series = pd.Series(
        uniprot_ensembl_df['To'].unique())
    uniprot_ensembl_unique_series.to_csv(os.path.join(
        out_dir, 'Uniprot_ensembl_mapping.txt'), index=False, header=False)


if __name__ == "__main__":
    args = parse_arguments()

    uniprot_ensembl_tsv = args.uniprot_ensembl_tsv

    main(uniprot_ensembl_tsv)
