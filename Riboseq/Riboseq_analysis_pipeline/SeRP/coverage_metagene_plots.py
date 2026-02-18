import os
import os.path

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO

import seaborn as sbn
from scipy.ndimage import uniform_filter1d


reference_fasta = '/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta'
SeRP_trans_dedup_dir = '/projects/splitorfs/work/Riboseq/Output/SeRP/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/deduplicated/filtered'


transcript_of_interest = ['ENST00000258962.5']
transcript_len_dict = {}
# obtain the length of the transcripts
for record in SeqIO.parse(reference_fasta, "fasta"):
    if record.id in transcript_of_interest:
        transcript_len_dict[record.id] = len(record.seq)


bam_files = []

for filename in os.listdir(SeRP_trans_dedup_dir):
    if filename.endswith("dedup_q10_filtered.bam"):
        bam_files.append(os.path.join(SeRP_trans_dedup_dir, filename))


Input_bam_files_CHX = [bam for bam in bam_files if 'In_CHX' in bam]
Input_bam_files_Puro = [bam for bam in bam_files if 'In_Puro' in bam]

IP_bam_files_CHX = [bam for bam in bam_files if 'IP_CHX' in bam]
IP_bam_files_Puro = [bam for bam in bam_files if 'IP_Puro' in bam]


for transcript in transcript_len_dict.keys():
    signal_track = np.zeros(
        (len(Input_bam_files_CHX) + 1, transcript_len_dict[transcript]), dtype=int)
    signal_track[0, :] = range(0, transcript_len_dict[transcript])
    for i, bam_file in enumerate(Input_bam_files_CHX):
        infile = pysam.AlignmentFile(bam_file, "rb")

        for alignment in infile.fetch(transcript):
            signal_track[i + 1, alignment.reference_start] += 1

    x = signal_track[0, :]
    df = pd.DataFrame({
        # repeat x for each y line
        "position on transcript in bp": np.tile(x, 3),
        "RPFs mapping": np.concatenate(signal_track[1:4, :]),
        "replicate": np.repeat(["In_CHX_1", "In_CHX_2", "In_CHX_4"], len(x))
    })

    sbn.lineplot(data=df, x="position on transcript in bp",
                 y="RPFs mapping", orient="y")
    sbn.lineplot(data=df, x="position on transcript in bp",
                 y="RPFs mapping", hue="replicate")

    signal_track_IP = np.zeros(
        (len(IP_bam_files_CHX) + 1, transcript_len_dict[transcript]), dtype=int)
    signal_track_IP[0, :] = range(0, transcript_len_dict[transcript])
    for i, bam_file in enumerate(IP_bam_files_CHX):
        infile = pysam.AlignmentFile(bam_file, "rb")

        for alignment in infile.fetch(transcript):
            signal_track_IP[i + 1, alignment.reference_start] += 1

    x = signal_track_IP[0, :]
    df_IP = pd.DataFrame({
        # repeat x for each y line
        "position on transcript in bp": np.tile(x, 3),
        "RPFs mapping": np.concatenate(signal_track_IP[1:4, :]),
        "replicate": np.repeat(["IP_CHX_1", "IP_CHX_2", "IP_CHX_4"], len(x))
    })

    sbn.lineplot(data=df_IP, x="position on transcript in bp",
                 y="RPFs mapping", orient="y")
    sbn.lineplot(data=df_IP, x="position on transcript in bp",
                 y="RPFs mapping", hue="replicate")

    signal_track_IP = signal_track_IP.astype(float)
    signal_track = signal_track.astype(float)

    signal_track_IP[1:4, :] = uniform_filter1d(
        signal_track_IP[1:4, :], size=25)
    signal_track[1:4, :] = uniform_filter1d(signal_track[1:4, :], size=25)
    signal_track_IP[1:4, :] += 1
    signal_track[1:4, :] += 1

signal_track_IP_over_In = np.log(
    np.divide(signal_track_IP[1:4, :], signal_track[1:4, :]))
df_IP_over_In_CHX = pd.DataFrame({
    # repeat x for each y line
    "position on transcript in bp": np.tile(x, 3),
    "log RPFs IP/In": np.concatenate(signal_track_IP_over_In),
    "replicate": np.repeat(["IP_over_In_CHX_1", "IP_over_In_CHX_2", "IP_over_In_CHX_4"], len(x))
})

sbn.lineplot(data=df_IP_over_In_CHX, x="position on transcript in bp",
             y="log RPFs IP/In", orient="y")

sbn.lineplot(data=df_IP_over_In_CHX, x="position on transcript in bp",
             y="log RPFs IP/In", hue="replicate")
