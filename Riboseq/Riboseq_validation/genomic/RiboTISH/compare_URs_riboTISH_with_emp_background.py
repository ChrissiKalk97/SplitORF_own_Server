import sys
import os
import pandas as pd
from Bio import SeqIO
import glob
import re

RiboTISH_out_dir = sys.argv[1]

# os.chdir("/projects/splitorfs/work/Riboseq/Output/RiboTISH_NMD_custom")
os.chdir(RiboTISH_out_dir)


for pred_ORFs_file in glob.glob("*.csv"):
    sample_name = pred_ORFs_file.rsplit('_', maxsplit = 3)[0]
    print(sample_name)
    RiboTISH_SO_df = pd.read_csv(pred_ORFs_file, header = 0)
    # print(RiboTISH_SO_df.head())
    RiboTISH_SO_df['ID'] = RiboTISH_SO_df['Gid'] + '|' + RiboTISH_SO_df['Tid']
    RiboTISH_SO_df['start'] = RiboTISH_SO_df['GenomePos'].apply(lambda x: int(re.split(r'[\:\-]', x)[1]))
    RiboTISH_SO_df['stop'] = RiboTISH_SO_df['GenomePos'].apply(lambda x: int(re.split(r'[\:\-]', x)[2]))
    # print(RiboTISH_SO_df[['start', 'stop']])
    SO_empirical_df = pd.read_csv(f'/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample/NMD_genome/{sample_name}_unique_regions.csv',
                                   header = 0)
    # SO_empirical_df = SO_empirical_df[SO_empirical_df['significant'] == 1]
    SO_empirical_df['ID'] = SO_empirical_df['name'].apply(lambda x: x.split(':')[0])
    SO_empirical_df = SO_empirical_df[['chr_unique', 'start', 'stop', 'name', 'ID', 'relative_count', 'num_reads', 'significant']]
    RiboTISH_SO_IDs = RiboTISH_SO_df['ID'].to_list()
    SO_empirical_df_filtered = SO_empirical_df[SO_empirical_df['ID'].isin(RiboTISH_SO_IDs)]
    # check data type of start and stop
    merged_df = pd.merge(SO_empirical_df_filtered, RiboTISH_SO_df, on='ID', how='outer')
    merged_df['same_region_start'] = merged_df.apply(lambda x: abs(x['start_y'] - x['start_x']) <= 3,axis = 1)
    merged_df['same_region_stop'] = merged_df.apply(lambda x: abs(x['stop_y'] - x['stop_x']) <= 3,axis = 1)
    merged_df['UR_within_mod_region'] = merged_df.apply(lambda x: (x['stop_y'] >= x['stop_x']) and (x['start_y'] <= x['start_x']),axis = 1)
    # filter out unique regions in the same transcript that were not found with Riboseq data
    merged_df = merged_df[merged_df['UR_within_mod_region'] == True]
    # with this filter we do not loose any hits
    merged_df = merged_df[merged_df['InFrameCount'] > 2]
    print('Nr URs validated:', len(merged_df.index))
    print('Nr unique URs:', merged_df['GenomePos'].nunique())
    merged_df_significant = merged_df[merged_df['significant'] == 1]
    print('Nr significant with empiricial approach:', merged_df_significant['GenomePos'].nunique())
    merged_df_insignificant = merged_df[merged_df['significant'] == 0]
    print('Nr non-significant with empiricial approach:', merged_df_insignificant['GenomePos'].nunique())

