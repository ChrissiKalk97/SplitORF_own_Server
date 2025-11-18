import sys
import pandas as pd
from pygtftk.gtf_interface import GTF
from Bio import SeqIO

# chae_gtf_path = sys.argv[1]
# transcript_fasta = sys.argv[2]
# out_file_name = sys.argv[3]

transcript_fasta = '/projects/serp/work/references/Chaetomium_thermophilum_transcript.fasta'
chae_gtf_path = '/projects/serp/work/references/Supplementary_File_2.gtf'
out_name = '/projects/serp/work/references/Chaetomium_thermophilum_'

chae_gtf_object = GTF(
    chae_gtf_path, check_ensembl_format=False)

# chae_gtf_object_rrna = chae_gtf_object.select_by_key('feature', 'transcript')\
#     .select_by_key('Old_gene_biotype', 'rRNA')

# chae_gtf_object_trna = chae_gtf_object.select_by_key('feature', 'transcript')\
#     .select_by_key('Old_gene_biotype', 'tRNA')

chae_gtf_object_protcod = chae_gtf_object.select_by_key('feature', 'transcript')\
    .select_by_key('gene_biotype', 'coding')

chae_gtf_object_noncod = chae_gtf_object.select_by_key('feature', 'transcript')\
    .select_by_key('gene_biotype', 'noncoding')

# rrna_transcripts = chae_gtf_object_rrna.get_tx_ids(nr=True)
# trna_transcripts = chae_gtf_object_trna.get_tx_ids(nr=True)
noncod_transcripts = chae_gtf_object_noncod.get_tx_ids(nr=True)
protcod_transcripts = chae_gtf_object_protcod.get_tx_ids(nr=True)

# with open(out_name + 'rRNA.fasta', "w") as out_handle:
#     # Parse input FASTA file
#     for record in SeqIO.parse(transcript_fasta, "fasta"):
#         # Check if any keyword is in the header
#         if any(word in record.description for word in rrna_transcripts):
#             SeqIO.write(record, out_handle, "fasta")


# with open(out_name + 'tRNA.fasta', "w") as out_handle:
#     # Parse input FASTA file
#     for record in SeqIO.parse(transcript_fasta, "fasta"):
#         # Check if any keyword is in the header
#         if any(word in record.description for word in trna_transcripts):
#             SeqIO.write(record, out_handle, "fasta")

with open(out_name + 'noncoding.fasta', "w") as out_handle:
    # Parse input FASTA file
    for record in SeqIO.parse(transcript_fasta, "fasta"):
        # Check if any keyword is in the header
        if record.description in noncod_transcripts:
            SeqIO.write(record, out_handle, "fasta")


with open(out_name + 'protein_coding.fasta', "w") as out_handle:
    # Parse input FASTA file
    for record in SeqIO.parse(transcript_fasta, "fasta"):
        # Check if any keyword is in the header
        if record.description in protcod_transcripts:
            SeqIO.write(record, out_handle, "fasta")
