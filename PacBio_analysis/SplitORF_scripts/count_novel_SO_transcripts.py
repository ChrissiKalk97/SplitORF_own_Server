import os
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    parser.add_argument("unique_protein_pairs_txt",
                        help="Path to UniqueProteinPairs.txt")

    parser.add_argument("novel_isoforms_txt",
                        help="Path to novel_isoforms.txt")

    parser.add_argument("--assembly_type",
                        help="assembly type: filtered or full")

    return parser.parse_args()


unique_protein_pairs_txt = '/home/ckalk/tools/SplitORF_pipeline/Output/run_12.09.2025-17.51.04_HUVEC_tama_merged/UniqueProteinORFPairs.txt'
novel_isoforms_txt = '/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama/HUVEC/compare_Ens_full_ref/HUVEC_merged_tama_gene_id_novel_isoforms.txt'


def main(unique_protein_pairs_txt, novel_isoforms_txt, assembly_type):

    outdir = os.path.dirname(novel_isoforms_txt)
    outname = os.path.basename(novel_isoforms_txt)[:-4]

    unique_protein_pairs_df = pd.read_csv(unique_protein_pairs_txt, sep='\t')
    split_orf_transcript_list = unique_protein_pairs_df['OrfTransID'].to_list()

    assert len(split_orf_transcript_list) == len(
        unique_protein_pairs_df['OrfTransID'].unique())

    novel_isoforms_df = pd.read_csv(novel_isoforms_txt)

    novel_split_orf_transcripts = novel_isoforms_df[novel_isoforms_df['qry_id'].isin(
        split_orf_transcript_list)]

    novel_split_orf_transcripts.to_csv(
        os.path.join(outdir, f'{outname}_splitorf.txt'))

    print(f'Number of novel Split-ORF transcripts compared to {assembly_type} assembly: ', len(
        novel_split_orf_transcripts.index))


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    unique_protein_pairs_txt = args.unique_protein_pairs_txt
    novel_isoforms_txt = args.novel_isoforms_txt
    assembly_type = args.assembly_type

    main(unique_protein_pairs_txt, novel_isoforms_txt, assembly_type)
