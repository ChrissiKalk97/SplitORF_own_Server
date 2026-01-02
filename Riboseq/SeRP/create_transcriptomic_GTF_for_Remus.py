# Create a pseudo transcriptomic GTF file for visulaization in IGV of transcriptomic
# mapped reads with CDS regions within transcripts

# usage: python

import os
import os.path
import argparse
import pandas as pd
import math


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create transcriptomic GTF.")

    parser.add_argument(
        "transcript_lengths_txt",
        type=str,
        help="Path to TXT file of transcript lengths"
    )

    parser.add_argument(
        "cds_coords_bed",
        type=str,
        help="Path to BED file of cds coordinates"
    )

    return parser.parse_args()


def main(transcript_lengths_txt, cds_coords_bed):
    transcript_lengths_df = pd.read_csv(transcript_lengths_txt, sep='\t')
    cds_coords_df = pd.read_csv(cds_coords_bed, sep='\t', header=None, names=[
                                'Transcript stable ID version', 'CDS_start', 'CDS_end'])

    cds_coords_df['Transcript stable ID'] = cds_coords_df['Transcript stable ID version'].apply(
        lambda x: x.split('.')[0])

    transcript_version_dict = dict(zip(
        cds_coords_df['Transcript stable ID'], cds_coords_df['Transcript stable ID version']))
    transcript_lengths_df['Transcript stable ID version'] = transcript_lengths_df['Transcript stable ID'].map(
        transcript_version_dict)

    cds_coords_df = cds_coords_df.set_index('Transcript stable ID version')

    gtf_string = ''

    for index, transcript in transcript_lengths_df.iterrows():
        transcript_length = transcript['Transcript length (including UTRs and CDS)']
        transcript_id = transcript['Transcript stable ID version']
        gene_id = transcript['Gene stable ID']

        if transcript_id in cds_coords_df.index:
            gtf_string += f'{transcript_id}\thavana\ttranscript\t1\t{transcript_length}\t.\t+\t.\tgene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_type "protein_coding"; transcript_type "protein_coding";\n'
            gtf_string += f'{transcript_id}\thavana\texon\t1\t{transcript_length}\t.\t+\t.\tgene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;\n'

            # am using the BED file so do need to add 1 to start
            cds_start = int(cds_coords_df.loc[transcript_id, 'CDS_start']) + 1
            cds_end = cds_coords_df.loc[transcript_id, 'CDS_end']
            if not math.isnan(float(cds_start)):
                gtf_string += f'{transcript_id}\thavana\tCDS\t{str(cds_start)}\t{cds_end}\t.\t+\t.\tgene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;\n'
            else:
                print('NA start')
                print(transcript_id)
                print(cds_coords_df.loc[transcript_id, 'CDS_start'])
                print(cds_start)
        else:
            print('Transcript not in CDS coord file')
            print(transcript['Transcript stable ID'])

    outdir = os.path.dirname(transcript_lengths_txt)
    transcriptomic_gtf = f'{outdir}/MANE_transcriptomic_GTF_for_Remus.gtf'

    with open(transcriptomic_gtf, "w") as gtf:
        gtf.write(gtf_string)


if __name__ == "__main__":
    args = parse_arguments()

    transcript_lengths_txt = args.transcript_lengths_txt
    cds_coords_bed = args.cds_coords_bed

    # cds_coords_bed = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/bowtie1/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed'
    # transcript_lengths_txt = '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/bowtie1/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_transcript_length.txt'

    main(transcript_lengths_txt, cds_coords_bed)
