import sys
import os
import pandas as pd
from Bio import SeqIO
import glob

os.chdir("/projects/splitorfs/work/Riboseq/Output/RiboTISH_NMD_custom")


for pred_ORFs_file in glob.glob("*_framebest_all.txt"):
    sample_name = pred_ORFs_file.rsplit('_', maxsplit = 5)[0]
    print(sample_name)
    SO_df= pd.read_csv(pred_ORFs_file, sep = "\t", header = 0)
    print('Nr of unique URs predicted', SO_df['GenomePos'].nunique())
    print(SO_df.sort_values('RiboPvalue')[['Gid', 'Tid', 'Symbol', 'GenomePos', 'TisType', 'RiboPvalue']])
    SO_df = SO_df.groupby('GenomePos').agg('first')
    SO_df = SO_df.reset_index()
    SO_df.sort_values('RiboPvalue')[['Gid', 'Tid', 'Symbol', 'GenomePos', 'TisType', 'RiboPvalue']].to_csv(f'{sample_name}_UR_translation_RiboTISH.csv', index=False)