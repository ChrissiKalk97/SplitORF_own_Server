# ------------------ IMPORTS ------------------ #

import sys

import pandas as pd
import argparse
import os
import os.path


sys.path.append(os.path.abspath("/home/ckalk/scripts/FuFis/src")) # noqa: F401

import GTF_Processing # noqa: F401




def parse_args():
    parser = argparse.ArgumentParser(
        description="Check for disease genes in validated Split-ORFs of respective disease Ribo-seq."
    )

    # Required positional arguments
    parser.add_argument("--validated_so_csv",
                        help="Path to validated Split-ORF CSV file")
    parser.add_argument("--disgenet_excel",
                        help="Path to disgenet Excel of respective disease")
    parser.add_argument("--samples",
                        help="Comma-sep string of samples to consider")
    parser.add_argument("--result_dir",
                        help="Path to directory for results")
    parser.add_argument("--region_type",
                        help="Region type, i.e. NMD or RI")
    parser.add_argument("--gtf",
                        help="GTF for Ensembl ID - Gene Symbol lifting")

    return parser.parse_args()


def main(validated_so_csv, disgenet_excel, region_type, result_dir, samples, gtf):
    # ------------------ DATA IMPORT ------------------ #
    disease = os.path.basename(os.path.dirname(disgenet_excel))

    validated_so_df = pd.read_csv(validated_so_csv)
    disease_df = pd.read_excel(disgenet_excel)
    gene_ids = list(validated_so_df['geneID'].unique())
    mapped_ids_dict, missed_ids = GTF_Processing.match_gene_identifiers(gene_ids, gtf_file=gtf, species='human',
                                                                        fields="symbol")

    validated_gene_names_df = {id: mapping['symbol']
                               for id, mapping in mapped_ids_dict.items()}
    validated_so_df['geneSymbol'] = validated_so_df['geneID'].map(
        validated_gene_names_df)

    sample_columns = [
        col for col in validated_so_df.columns if col.startswith(samples)]
    columns_keep = sample_columns + ['geneSymbol', 'geneID', 'OrfTransID',
                                     'OrfIDs', 'OrfStarts', 'OrfStart', 'OrfIndex', 'OrfPosition']

    validated_so_df = validated_so_df[columns_keep]
    validated_so_df[f'ValidationCount_{samples}'] = validated_so_df.iloc[:, 0:len(
        sample_columns)].apply(lambda x: sum(x), axis=1)

    validated_so_df = validated_so_df[validated_so_df[f'ValidationCount_{samples}'] > 0]
    
    validated_so_df.to_csv(os.path.join(
        result_dir, f'{region_type}_{samples}_validated_regions.csv'))
    
    print(f'Number of validated genes for {samples}:', len(validated_so_df.index))

    validated_so_df = validated_so_df[validated_so_df['geneSymbol'].isin(
        disease_df['gene_symbol'])]

    disease_df = disease_df[disease_df['gene_symbol'].isin(
        validated_so_df['geneSymbol'])]

    validated_so_df.to_csv(os.path.join(
        result_dir, f'{region_type}_{samples}_{disease}_DISGENET_intersection.csv'))
    
    disease_df.to_csv(os.path.join(
        result_dir, f'{disease}_{region_type}_{samples}_Riboseq_intersection.csv'))


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    validated_so_csv = args.validated_so_csv
    disgenet_excel = args.disgenet_excel
    samples = args.samples
    result_dir = args.result_dir
    region_type = args.region_type
    gtf = args.gtf

    # validated_so_csv = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/validated_so_df.csv'
    # disgenet_excel = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/DISGENET_analysis/glioblastoma/search_result_C0017636-C1621958-C1514422-C0278878-C0349543+5.xlsx'
    # samples = 'SRR10'
    # region_type = 'NMD'
    # gtf = '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf'

    main(validated_so_csv, disgenet_excel,
         region_type, result_dir, samples, gtf)
