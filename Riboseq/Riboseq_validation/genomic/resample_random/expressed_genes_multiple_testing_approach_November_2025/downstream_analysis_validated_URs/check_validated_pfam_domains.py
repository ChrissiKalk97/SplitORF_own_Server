# This script takes the validated unique regions and checks for their annotated
# PFAM domains
# ------------------ IMPORTS ------------------ #
import os
import os.path
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description='Summarize the trans.'
    )

    # Required positional arguments

    parser.add_argument('--annotated_pfam_file',
                        help='Path to the Split-ORF pipeline output with annotated PFAM domains')
    parser.add_argument('--validation_file',
                        help='of the validated_so_df.csv which is checked for RBP validations')
    parser.add_argument('--region_type',
                        help='NMD or RI')
    parser.add_argument('--rbp_file',
                        help='Path to RNA-binding protein CSV from RBPDB')
    parser.add_argument('--pfam_annotation',
                        help='Path to PFAM annotation file: Pfam-A.clans.tsv')

    return parser.parse_args()


def main(annotated_pfam_file, validation_file, region_type, rbp_file, pfam_annotation):
    outdir = os.path.dirname(os.path.dirname(os.path.dirname(validation_file)))
    os.makedirs(os.path.join(outdir, 'PFAM_analysis'), exist_ok=True)
    os.makedirs(os.path.join(outdir, 'PFAM_analysis',
                region_type), exist_ok=True)

    so_annotated_pfam_df = pd.read_csv(annotated_pfam_file, sep='\t')

    pfam_annotation_df = pd.read_csv(pfam_annotation, sep='\t', header=None, names=[
                                     'PfamID', 'ClanID', 'ClanName', 'PfamName', 'PfamDescription'])

    validation_df = pd.read_csv(validation_file)
    validation_df = validation_df.iloc[:, 1:].copy()

    # which validated Split-ORF transcripts have 2 or more PFAM domains?
    so_val_pfam_df = so_annotated_pfam_df[so_annotated_pfam_df['OrfTransID'].isin(
        validation_df['OrfTransID'])].copy()
    so_val_pfam_df['ORF-DomainAnnot_list'] = so_val_pfam_df['ORF-DomainAnnot'].apply(
        lambda x: str(x).split(',')).copy()

    # add clan information to the SO df with PFAM info
    so_val_pfam_df_exploded = so_val_pfam_df.explode('ORF-DomainAnnot_list')
    so_val_pfam_df_exploded_clan_info = so_val_pfam_df_exploded.merge(
        pfam_annotation_df, left_on='ORF-DomainAnnot_list', right_on='PfamID', how='left')
    list_cols = pfam_annotation_df.columns.tolist() + ['ORF-DomainAnnot_list']
    agg_dict = {col: (list if col in list_cols else "first")
                for col in so_val_pfam_df_exploded_clan_info.columns if col != 'OrfTransID'}
    so_val_pfam_df = so_val_pfam_df_exploded_clan_info.groupby(
        'OrfTransID').agg(agg_dict).reset_index()

    # subset for valdiated SOs with 1 PFAM domain or more
    so_val_with_pfam_df = so_val_pfam_df[so_val_pfam_df['ORF-DomainAnnot'].notna()]

    # subset for valdiated SOs with 2 PFAM domains or more
    so_val_with_two_pfams_df = so_val_with_pfam_df[so_val_with_pfam_df['ORF-DomainAnnot_list'].apply(
        lambda x: len(x) > 1)]

    print(f'Number of Split-ORF genes validated {region_type}:',
          len(so_val_pfam_df['geneID'].unique()))
    print(f'Number of Split-ORF genes validated with PFAM domains {region_type}:',
          len(so_val_with_pfam_df['geneID'].unique()))
    print(f'Number of Split-ORF genes validated with 2 or more PFAM domains {region_type}:',
          len(so_val_with_two_pfams_df['geneID'].unique()))

    # write results SOs with PFAMs
    so_val_with_pfam_df.to_csv(os.path.join(outdir, 'PFAM_analysis',
                                            region_type, 'so_val_with_pfam_df.csv'))
    (
        so_val_with_pfam_df['ORF-DomainAnnot_list']
        .explode()
        .dropna()
        .drop_duplicates()
        .to_csv(os.path.join(outdir, 'PFAM_analysis',
                             region_type, 'validated_pfams.csv'), index=False, header=False)
    )
    # clan analysis for SOs with PFAMs, count each PFAM domain once per gene
    (
        so_val_with_pfam_df[['ClanID', 'geneID']]
        .explode(['ClanID'])
        .drop_duplicates()['ClanID']
        .explode()
        .value_counts()
        .to_csv(os.path.join(outdir, 'PFAM_analysis', region_type, 'validated_pfam_clans.csv'))
    )

    # write results SOs with 2 PFAMs
    so_val_with_two_pfams_df.to_csv(os.path.join(outdir, 'PFAM_analysis',
                                                 region_type, 'so_val_with_two_pfams_df.csv'))
    (
        so_val_with_two_pfams_df['ORF-DomainAnnot_list']
        .explode()
        .dropna()
        .drop_duplicates()
        .to_csv(os.path.join(outdir, 'PFAM_analysis',
                             region_type, 'validated_pfams_of_trans_with_2_or_more.csv'), index=False, header=False)
    )

    # Which val transcripts have PFAM domains and are RBPs as well?
    rbp_df = pd.read_csv(rbp_file, sep=';')
    rbp_df = rbp_df.iloc[:, 1:].copy()
    rbp_df_filtered = rbp_df[rbp_df['Annotation ID'].isin(
        validation_df['geneID'])].copy()
    rbp_df_filtered = rbp_df_filtered.rename(
        columns={'Annotation ID': 'geneID'})

    rbp_val_with_pfam_df = so_val_with_pfam_df.merge(
        rbp_df_filtered, on='geneID', how='inner')

    rbp_val_with_pfam_df.to_csv(os.path.join(outdir, 'PFAM_analysis',
                                             region_type, 'rbp_val_with_pfam_df.csv'))

    (
        rbp_val_with_pfam_df['ORF-DomainAnnot_list']
        .explode()
        .dropna()
        .drop_duplicates()
        .to_csv(os.path.join(outdir, 'PFAM_analysis',
                             region_type, 'validated_rbp_pfams_of_trans_with_2_or_more.csv'), index=False, header=False)
    )
    # clan analysis for SOs with PFAMs
    # clan analysis for SOs with PFAMs, count each PFAM domain once per gene
    (
        rbp_val_with_pfam_df[['ClanID', 'geneID']]
        .explode(['ClanID'])
        .drop_duplicates()['ClanID']
        .explode()
        .value_counts()
        .to_csv(os.path.join(outdir, 'PFAM_analysis', region_type, 'validated_rbp_pfam_clans.csv'))
    )

    print(f'Number of RBPs with PFAM domain {region_type}:', len(
        rbp_val_with_pfam_df['geneID'].unique()))

    # Which val transcripts have 2 PFAM domains and are RBPs as well?
    rbp_val_with_two_pfams_df = so_val_with_two_pfams_df.merge(
        rbp_df_filtered, on='geneID', how='inner')

    rbp_val_with_two_pfams_df.to_csv(os.path.join(outdir, 'PFAM_analysis',
                                                  region_type, 'rbp_val_with_two_pfams_df.csv'))

    (
        rbp_val_with_two_pfams_df['ORF-DomainAnnot_list']
        .explode()
        .dropna()
        .drop_duplicates()
        .to_csv(os.path.join(outdir, 'PFAM_analysis',
                             region_type, 'validated_rbp_pfams.csv'), index=False, header=False)
    )

    print(f'Number of RBPs with 2 or more PFAM domains {region_type}:', len(
        rbp_val_with_two_pfams_df['geneID'].unique()))


if __name__ == '__main__':
    args = parse_args()

    annotated_pfam_file = args.annotated_pfam_file
    validation_file = args.validation_file
    rbp_file = args.rbp_file
    region_type = args.region_type
    pfam_annotation = args.pfam_annotation

    # annotated_pfam_file = '/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-13.07.17_NMD_for_paper/UniqueProteinORFPairs_annotated.txt'
    # validation_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/validated_so_df.csv'
    # region_type = 'NMD'
    # rbp_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/proteins.php'
    # pfam_annotation = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/PFAM_analysis/Pfam-A.clans.tsv'

    main(annotated_pfam_file, validation_file,
         region_type, rbp_file, pfam_annotation)
