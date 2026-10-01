#######################################################################
## get the protein coding and NMD transcripts from fiftynt rule and  ##
## ORFanage FIRST mode                                               ##
#######################################################################
# usage

# python get_nmd_protein_coding_transcripts.py Orfanage_50nt.csv custom_50nt.csv cell_type output_dir

import os
import sys
from typing import List
import pandas as pd

# ORFanage_CDS_file = '~/tools/NMD_fetaure_composition/Output/CM_TAMA_ORFanage_FIRST_09_09_26/CM_TAMA_ORFanage_FIRST_09_09_26.csv'
# custom_CDS_file = '~/tools/NMD_fetaure_composition/Output/CM_merged_tama_10000_10000_iso_mando_stringtie_50nt/CM_merged_tama_10000_10000_iso_mando_stringtie_50nt.csv'
# ORFanage_CDS = pd.read_csv(ORFanage_CDS_file, header=0, index_col=0)
# # custom CDS prediction
# custom_CDS = pd.read_csv(custom_CDS_file, header=0, index_col=0)

# csv file of the ORFanage CDS prediction (for ease of data access, start ORF, end ORF)
ORFanage_CDS = pd.read_csv(sys.argv[1], header=0, index_col=0)
# custom CDS prediction
custom_CDS = pd.read_csv(sys.argv[2], header=0, index_col=0)

cell_type = sys.argv[3]
output_dir = sys.argv[4]


orfanage_50nt = set(ORFanage_CDS[ORFanage_CDS['50_nt'] == 1].index)
custom_50nt = set(custom_CDS[custom_CDS['50_nt'] == 1].index)

all_50nt_transcripts = orfanage_50nt | custom_50nt

orfanage_protein_coding = set(ORFanage_CDS[(ORFanage_CDS['50_nt'] == 0) & (
    ~ORFanage_CDS.index.isin(all_50nt_transcripts))].index)

with open(os.path.join(output_dir, f'{cell_type}_protein_coding_ORFanage_FIRST.txt'), 'w') as f:
    for transcript in orfanage_protein_coding:
        f.write(f"{transcript}\n")


with open(os.path.join(output_dir, f'{cell_type}_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline.txt'), 'w') as f:
    for transcript in all_50nt_transcripts:
        f.write(f"{transcript}\n")
