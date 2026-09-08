# create BED file from the Excel of publicaiton

# conda: SplitORF

import sys
import os.path
import pandas as pd


uce_excel = "/projects/splitorfs/work/split-orf-prediction/ultraconserved_regions/42003_2025_9115_MOESM4_ESM.xlsx"


uce_df = pd.read_excel(uce_excel, sheet_name=4)

exonic_uce_df = uce_df.iloc[:, 0:4].dropna().copy()
intronic_uce_df = uce_df.iloc[:, 5:].dropna().copy()


exonic_uce_df.columns = exonic_uce_df.iloc[0]  # use first for as header
exonic_uce_df = exonic_uce_df[1:].reset_index(drop=True)   # drop row


intronic_uce_df.columns = intronic_uce_df.iloc[0]  # use first for as header
intronic_uce_df = intronic_uce_df[1:].reset_index(drop=True)   # drop row

intronic_uce_df['name'] = intronic_uce_df['name'] + '_intronic'
exonic_uce_df['name'] = exonic_uce_df['name'] + '_exonic'

merged_uce_df = pd.concat([exonic_uce_df, intronic_uce_df], ignore_index=True)

merged_uce_df[["chr", "start", "end", "name"]].to_csv(
    "/projects/splitorfs/work/split-orf-prediction/ultraconserved_regions/UCE_regions.bed", sep="\t", header=False, index=False
)
