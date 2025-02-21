"""
Script Name: prepare_unique_positions_as_ORF_input.py
Description: ...
.
Usage: python prepare_unique_positions_as_ORF_input.py unique_transcript_regions.bed output_dir
unique_transcript_regions.bed = '/projects/splitorfs/work/Riboseq/data/region_input/transcriptomic/Unique_DNA_Regions_for_riboseq_NMD.bed'
output_dir = '/projects/splitorfs/work/Riboseq/data/RiboTISH'
"""



import sys
import os
import pandas as pd
from pybedtools import BedTool
import os.path

unique_trans_region = sys.argv[1]
output_dir = sys.argv[2]


unique_trans_region_df = pd.read_csv(unique_trans_region, sep = '\t', header = None)
unique_trans_region_df.columns = ['ID', 'start', 'stop', 'ORF']

sum((unique_trans_region_df['stop'] - unique_trans_region_df['start']) % 3 == 0)
# 1633 out of 4530 do not need any polishing, the rest does

# 1. get ORF coords
# compare the start sites
# 3. if the starts are the same, add  3- what is left from modulo
# if they differ: remove from the start 3 minus what is left from the modulo
# IDEA: extend the unique regions if not a multiple of 3, this should of course not be performed 
# outside of tHe ORF

# calculate the modulo as a new column
unique_trans_region_df['modulo'] = (unique_trans_region_df['stop'] - unique_trans_region_df['start']) % 3
# get the ORF start as a new column
unique_trans_region_df['ORF_start'] = unique_trans_region_df['ORF'].apply(lambda x: x.split(':')[1])

def create_region_for_RiboTISH(ORF_start, region_start, region_end):
    ORF_start = int(ORF_start)
    start_modulo = (region_start - ORF_start) % 3
    end_modulo = (region_end - ORF_start) % 3
    if start_modulo == 0:
        return [region_start, region_end + (3 - end_modulo)]
    else:
        if end_modulo > 0:
            return [region_start - start_modulo, region_end + (3 - end_modulo)]
        else: 
            return [region_start - start_modulo, region_end]

unique_trans_region_df['modified_regions'] = unique_trans_region_df.apply(lambda x: create_region_for_RiboTISH(x['ORF_start'],  x['start'], x['stop']), axis = 1)


def create_BedTool_string(ID, modified_regions):
    BedTool_string = ''
    BedTool_string = BedTool_string + ID.split('|')[1] + '\t'
    BedTool_string = BedTool_string + str(modified_regions[0]) + '\t' + str(modified_regions[1]) + '\n'
    return BedTool_string

unique_trans_region_df['bedtool_string'] = unique_trans_region_df.apply(lambda x: create_BedTool_string(x['ID'], x['modified_regions']), axis = 1)

bed_string = ''.join(unique_trans_region_df['bedtool_string'].to_list())
mod_unique_region_bedtool = BedTool(bed_string, from_string = True)
out_prefix = os.path.basename(unique_trans_region)
mod_unique_region_bedtool.saveas(output_dir + '/' + out_prefix + '_RiboTISH_modified.bed')