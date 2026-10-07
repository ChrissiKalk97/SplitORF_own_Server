#######################################################################
## analyze_gtf_vs_csv_ORF_predictions.py  compares the ORFanage CDS  ##
## predictions with the custom ones in the csv file from the feat    ##
## calculation pipeline                                              ##
#######################################################################
# usage

# python analyze_gtf_vs_csv_ORF_predictions.py ../Orfanage_long_trans_CDS.gtf /Users/christina/Documents/NMD_prediction/fiftynt_rule/Output/long_assembly_to_compare_to_orfanage/long_assembly_to_compare_to_orfanage.csv


import sys
from typing import List
import pandas as pd
from pygtftk.gtf_interface import GTF

# csv file of the ORFanage CDS prediction (for ease of data access, start ORF, end ORF)
ORFanage_CDS = pd.read_csv(sys.argv[1], header=0, index_col=0)
# custom CDS prediction
custom_CDS = pd.read_csv(sys.argv[2], header=0, index_col=0)

tids_custom = set(custom_CDS.index)
tids_ORFanage = set(ORFanage_CDS.index)
tids_intersection = set(tid for tid in tids_custom if tid in tids_ORFanage)
tids_only_ORF = set(tid for tid in tids_ORFanage if tid not in tids_custom)
tids_only_custom = set(tid for tid in tids_custom if tid not in tids_ORFanage)
tids_union = tids_custom | tids_ORFanage


print('Number of transcripts with ORF in union', len(tids_union))
print('Number of transcripts with ORF in ORFanage and custom',
      len(tids_intersection))
print('Number of transcripts with CDS only defined with custom', len(tids_only_custom))
print('Number of transcripts with CDS only defined with ORFanage', len(tids_only_ORF))
print('\n')

# How many of the custom unique transcripts are 50nt positive?
custom_only_50 = custom_CDS[(custom_CDS['50_nt'] == 1) & (
    custom_CDS.index.isin(tids_only_custom))]
print('The nr of 50nt positive transcripts for the unique custom transcripts is:',
      custom_only_50.shape[0])
print('\n')


# How many of the custom ORFs are 50nt positive and 0 or absent from ORFanage?
orfanage_50 = ORFanage_CDS[ORFanage_CDS['50_nt'] == 1].index
custom_CDS_50_no_orfanage = custom_CDS[(custom_CDS['50_nt'] == 1) & (
    ~custom_CDS.index.isin(orfanage_50))]
print('The nr of 50nt positive transcripts for the custom transcripts which are not present or 0 in ORFanage:',
      custom_CDS_50_no_orfanage.shape[0])
print('\n')


# How many of the ORfanage ORFs are 50nt positive but negative or absent from custom?
custom_50 = custom_CDS[custom_CDS['50_nt'] == 1].index
orfanage_CDS_50_no_custom = ORFanage_CDS[(ORFanage_CDS['50_nt'] == 1) & (
    ~ORFanage_CDS.index.isin(custom_50))]
print('The nr of 50nt positive transcripts for the ORfanage transcripts which are not present or 0 in custom:',
      orfanage_CDS_50_no_custom.shape[0])
print('\n')

# Step 1: Merge the two dataframes on their indices
merged_df = ORFanage_CDS.merge(
    custom_CDS, left_index=True, right_index=True, suffixes=('_ORF', '_custom'))

# Step 2: Compare the values in the specified column
comparison_column = 'end_ORF'
matches = merged_df[merged_df[f'{comparison_column}_ORF'].astype(
    int) == merged_df[f'{comparison_column}_custom'].astype(int)+4]
non_matches = merged_df[~(merged_df[f'{comparison_column}_ORF'].astype(
    int) == merged_df[f'{comparison_column}_custom'].astype(int)+4)]

# Step 3: Count the number of matches
number_of_matches = matches.shape[0]
print(
    f'The number of matching values in the column "{comparison_column}" for the common indices is: {number_of_matches}')
print(
    f'The percentage of matching values in the column "{comparison_column}" for the common indices is: {number_of_matches/len(merged_df.index)}')
print('\n')

# investigate how often the 50nt is the same for the match cases
comparison_column = '50_nt'
same_50 = matches[matches[f'{comparison_column}_ORF']
                  == matches[f'{comparison_column}_custom']]
ORF_50 = matches[(matches[f'{comparison_column}_ORF'] == 1) & (
    matches[f'{comparison_column}_custom'] == 0)]
custom_50 = matches[(matches[f'{comparison_column}_ORF'] == 0) & (
    matches[f'{comparison_column}_custom'] == 1)]

print(
    f'The number of times that the 50 nt status is same for same CDS end is {same_50.shape[0]}')
print(
    f'The number of times that custom had 50 plus and ORFanage 50 zero for same end is : {custom_50.shape[0]}')
print(
    f'The number of times that custom had 50 zero and ORFanage 50 plus for same end is : {ORF_50.shape[0]}')
print('\n')


# Investigate the non-matching end-coordinates: how often does each have 50nt plus, while the other one is minus?
non_ORF_50 = non_matches[(non_matches[f'{comparison_column}_ORF'] == 1) & (
    non_matches[f'{comparison_column}_custom'] == 0)]
non_custom_50 = non_matches[(non_matches[f'{comparison_column}_ORF'] == 0) & (
    non_matches[f'{comparison_column}_custom'] == 1)]
print(
    f'The number of times that custom had 50 plus and ORFanage 50 zero for DIFFERENT end is : {non_custom_50.shape[0]}')
print(
    f'The number of times that custom had 50 zero and ORFanage 50 plus for DIFFERENT end is : {non_ORF_50.shape[0]}')
print('\n')


# COMPARE START INDICES AS WELL
comparison_column = 'start_ORF'
matches = merged_df[merged_df[f'{comparison_column}_ORF'].astype(
    int) == merged_df[f'{comparison_column}_custom'].astype(int)+1]

# Step 3: Count the number of matches
number_of_matches = matches.shape[0]

print(
    f'The number of matching values in the column "{comparison_column}" for the common indices is: {number_of_matches}')
print(
    f'The percentage of matching values in the column "{comparison_column}" for the common indices is: {number_of_matches/len(merged_df.index)}')
