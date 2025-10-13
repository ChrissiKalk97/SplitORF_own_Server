import sys
import os
import os.path
import pandas as pd

DEG_csv = sys.argv[1]
tID_txt = sys.argv[2]
outname = sys.argv[3]

DEG_df = pd.read_csv(DEG_csv, index_col=0)
tID_Series = pd.read_csv(tID_txt, header=None)

DEG_df_filtered = DEG_df[DEG_df.index.isin(tID_Series[0].to_list())]

DEG_df_filtered.to_csv(outname)
