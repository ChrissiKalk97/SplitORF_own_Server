# Take any GTF file and extract the exon coordinates of the transcript
# in the format required by the Split-ORF pipeline

# usage get_exon_coords_from_gtf.py in.gtf out.bed


import argparse
import pandas as pd
from pygtftk.gtf_interface import GTF


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get exon coordinates of transcript for Split-ORF pipeline.")

    parser.add_argument(
        "gtf_file",
        type=str,
        help="Path to the input GTF file"
    )

    parser.add_argument(
        "output_bed_file",
        type=str,
        help="Path to the output file with the exon coords (bed)"
    )

    return parser.parse_args()


def main(path_to_custom_gtf, outfile):
    # path_to_custom_gtf = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/HUVEC/Orfanage/Orfanage_FIRST_09_09_26/HUVEC_TAMA_ORFanage_FIRST_CDS_numbered.gtf'
    # load custom gtf
    custom_gtf = GTF(path_to_custom_gtf, check_ensembl_format=False)

    # get transcript info
    mando_cds_info = custom_gtf.select_by_key(
        "feature", "CDS").extract_data("transcript_id,start,end,strand,chr,phase")

    mando_cds_info_df = mando_cds_info.as_data_frame()

    # get transcript info
    mando_transcript_info = custom_gtf.select_by_key(
        "feature", "transcript").extract_data("gene_id,transcript_id,start,end,strand,chr")

    mando_transcript_info_df = mando_transcript_info.as_data_frame()

    # merge dataframes based on transcript ID
    transcript_cds_df = pd.merge(mando_cds_info_df, mando_transcript_info_df,
                                 on='transcript_id', suffixes=('_cds', '_transcript'), how='outer')

    transcript_cds_df['start_cds'] = transcript_cds_df['start_cds'].astype(
        int)
    transcript_cds_df['start_transcript'] = transcript_cds_df['start_transcript'].astype(
        int)

    transcript_cds_df['gid|tid'] = transcript_cds_df['gene_id'] + \
        '|' + transcript_cds_df['transcript_id']

    # reorder columns to the correct order
    # note that each transcript in the GTF has a CDS!
    transcript_cds_df = transcript_cds_df[['seqid_transcript', 'start_cds',
                                           'end_cds', 'gid|tid', 'phase', 'strand_transcript']]

    transcript_cds_df.to_csv(outfile, sep='\t', index=False, header=False)


if __name__ == "__main__":
    args = parse_arguments()
    main(args.gtf_file, args.output_bed_file)
