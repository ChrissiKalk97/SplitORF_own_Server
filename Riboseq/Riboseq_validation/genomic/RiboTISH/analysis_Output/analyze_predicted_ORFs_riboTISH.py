import sys
import os
import pandas as pd
from Bio import SeqIO

os.chdir("/projects/splitorfs/work/Riboseq/Output/RiboTISH")
heart1 = pd.read_csv("ERR3367797_pred_framebest.txt", sep = "\t", header = 0)
# over 130 000 translating ORFs predicted in that sample

heart1.groupby("Tid").count()
# 107 558 different transcripts are predicted to harbour translating ORFs

# this is still sounds a bit much to me

heart1[heart1['RiboPvalue'] > 0.05]
# there are none that have a p-value above 5%, however, these are
# p-values not Q-vlaues, so they are not corrected

heart1[heart1['FrameQvalue'] > 0.05]
# there are no entries above 5%
# FDR filtering was coduncted

# read in NMD transcripts and extract transcript IDs
NMD_trans_dict = SeqIO.to_dict(SeqIO.parse("/projects/splitorfs/work/reference_files/Ensembl_transcripts/NMD_transcripts_CDNA.fa", "fasta"))
NMD_trans_IDs = list(NMD_trans_dict.keys())

heart1['Gid_Tid'] = heart1['Gid'] + '|' + heart1['Tid']
heart1_nmd = heart1[heart1['Gid_Tid'].isin(NMD_trans_IDs)]
# 24386 predicted translating ORFs
heart1_nmd.groupby('Gid_Tid').count()
# 15361 NMD transcripts predicted to be translating