import sys
from typing import List
import pandas as pd

# this was still the run with the incorrect Input (only ORFanage transcripts as protein coding transcripts...)
so_ens = '/projects/splitorfs/work/split-orf-prediction/Output/run_05.10.2026-11.57.15_HUVEC_assembly_Ens110_ref/UniqueProteinORFPairs.txt'
so_custom = '/projects/splitorfs/work/split-orf-prediction/Output/run_02.10.2026-09.50.31_HUVEC_assembly_ref_ORFanage_FIRST_and_Ens110/UniqueProteinORFPairs.txt'

so_ens_df = pd.read_csv(so_ens, sep='\t')
so_custom_df = pd.read_csv(so_custom, sep='\t')


so_ens = '/projects/splitorfs/work/split-orf-prediction/Output/run_06.10.2026-09.30.04_CM_assembly_Ens110_ref/UniqueProteinORFPairs.txt'
so_custom = '/projects/splitorfs/work/split-orf-prediction/Output/run_06.10.2026-12.47.25_CM_assembly_ref_ORFanage_FIRST_and_Ens110/UniqueProteinORFPairs.txt'

so_ens_df = pd.read_csv(so_ens, sep='\t')
so_custom_df = pd.read_csv(so_custom, sep='\t')


so_ens_ur = '/projects/splitorfs/work/split-orf-prediction/Output/run_06.10.2026-09.30.04_CM_assembly_Ens110_ref/Unique_DNA_Regions_transcriptomic.bed'
so_custom_ur = '/projects/splitorfs/work/split-orf-prediction/Output/run_06.10.2026-12.47.25_CM_assembly_ref_ORFanage_FIRST_and_Ens110/Unique_DNA_Regions_transcriptomic.bed'

so_ens_ur_df = pd.read_csv(so_ens_ur, sep='\t', header = None)
so_custom_ur_df = pd.read_csv(so_custom_ur, sep='\t', header = None)
so_ens_ur_df[~so_ens_ur_df.iloc[:,0].isin(so_custom_ur_df.iloc[:,0])]