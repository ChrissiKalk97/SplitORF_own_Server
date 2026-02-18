from pygtftk.gtf_interface import GTF
import argparse
from typing import List


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('in_gtf',
                        help='Path to Input GTF to filter')

    parser.add_argument('out_gtf',
                        help='Path to output GTF')

    return parser.parse_args()


def main(in_gtf, out_gtf):
    def get_transcript_string(transcript_ids: List[str]) -> str:
        '''get string of transcript or other ids for filtering 
        of a GTF class object from pygtftk'''
        transcript_string = ''
        if len(transcript_ids) > 1:
            for transcript_id in transcript_ids[:-1]:
                transcript_string += transcript_id+','
        transcript_string += transcript_ids[-1]
        return transcript_string

    # load GTF
    gtf_to_filter = GTF(in_gtf, check_ensembl_format=False)
    # get gene: transcript mapping
    gene_tx_dict = gtf_to_filter.get_gn_to_tx()
    # get transcript: length mapping
    transcript_length_dict = gtf_to_filter.get_transcript_size()

    # get longest transcript per gene
    longest_transcript_list = []
    for key, value in gene_tx_dict.items():
        if len(value) == 1:
            longest_transcript_list.append(value[0])
        else:
            gene_dict = {k: v for k,
                         v in transcript_length_dict.items() if k in value}
            longest_transcript_key = max(gene_dict, key=gene_dict.get)
            longest_transcript_list.append(longest_transcript_key)

    # check: are there as many transcripts as genes
    assert len(longest_transcript_list) == len(gene_tx_dict)

    # filter GTF for longest transcripts
    longest_transcripts_string = get_transcript_string(longest_transcript_list)
    longest_transcript_gtf = gtf_to_filter.select_by_key(
        'transcript_id', longest_transcripts_string)

    longest_transcript_gtf.write(out_gtf)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    in_gtf = args.in_gtf
    out_gtf = args.out_gtf

    # in_gtf = "/projects/serp/work/references/Supplementary_File_2.gtf"
    main(in_gtf, out_gtf)
