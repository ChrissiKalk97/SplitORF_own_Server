import sys
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

fasta_file = sys.argv[1]
gtf = sys.argv[2]
outfasta = sys.argv[3]

# fasta_file="/home/ckalk//tools/SplitORF_pipeline/Input2023/CM_assembly/CM_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID.fa"
# gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/CM/CM_LR_SR_support_filtered.gtf"
# outfasta="/projects/splitorfs/work/LLMs/TranslationAI/Input/CM_assembly/CM_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID_for_transAI.fa"

gtf_df = pd.read_csv(gtf, sep='\t', header=None, names=[
                     'chr', 'source', 'feature', 'start', 'stop', 'phase', 'strand', 'sth', 'information'], dtype={'chr': str})
gtf_transcript_df = gtf_df[gtf_df['feature'] == 'transcript'].copy()
gtf_transcript_df['gID'] = gtf_transcript_df['information'].apply(
    lambda x: x.split(';')[0].split('"')[1])
gtf_transcript_df['tID'] = gtf_transcript_df['information'].apply(
    lambda x: x.split(';')[1].split('"')[1])
gtf_transcript_df['gID|tID'] = gtf_transcript_df['gID'] + \
    '|' + gtf_transcript_df['tID']


with open(fasta_file) as input_handle, open(outfasta, 'w') as output_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        gtf_trans_filtered_df = gtf_transcript_df[gtf_transcript_df['gID|tID'] == record.id].copy(
        )
        record.id = gtf_trans_filtered_df['chr'].item() + ':' + str(gtf_trans_filtered_df['start'].item() - 1) + '-' + str(gtf_trans_filtered_df['stop'].item()) + \
            '(' + gtf_trans_filtered_df['strand'].item() + \
            ')' + '(' + record.id + ')' + '(10, 100,)'
        record.description = ''
        record.seq = Seq(str(record.seq).replace('\n', ''))
        seq_string = str(record.seq).replace('\n', '')
        output_handle.write(f">{record.id}\n{seq_string}\n")
