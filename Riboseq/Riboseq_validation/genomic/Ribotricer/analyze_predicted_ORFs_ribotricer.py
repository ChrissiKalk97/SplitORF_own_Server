import sys
import os
import pandas as pd

os.chdir("/projects/splitorfs/work/Riboseq/Output/Ribotricer")
OHMX20220060_001 = pd.read_csv("OHMX20220060_001_RI_sorted.bam_detected_translation_translating_ORFs.tsv", sep = "\t", header = 0)


OHMX20220060_001_filtered = OHMX20220060_001[OHMX20220060_001["valid_codons_ratio"] > 0.5]
OHMX20220060_001_filtered.groupby("transcript_id").count()

nmd_trans = OHMX20220060_001_filtered[OHMX20220060_001_filtered["transcript_type"] == "nonsense_mediated_decay"]
nmd_trans.groupby("transcript_id").count()
nmd_trans.groupby("transcript_id").count()["ORF_ID"].mean()
# this is still 3.35 for the NMD transcripts on average
# 2857 NMD transcripts left after the filtering with 0.5