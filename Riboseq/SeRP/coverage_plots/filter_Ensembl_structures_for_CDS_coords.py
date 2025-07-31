import sys
import pandas as pd
from pygtftk.gtf_interface import GTF

MANE_path = sys.argv[1]
MANE_GTF_path = sys.argv[2]
out_file_name = sys.argv[3]

Ensembl_MANE_114 = pd.read_csv(
    MANE_path, sep='\t')
MANE_GTF = GTF(
    MANE_GTF_path)


MANE_txids = MANE_GTF.get_tx_ids(nr=True)


Ensembl_Ids = list(Ensembl_MANE_114['Transcript stable ID version'].unique())

Ensembl_MANE_114 = Ensembl_MANE_114[Ensembl_MANE_114['Transcript stable ID version'].isin(
    MANE_txids)]

Ensembl_canonical_starts = Ensembl_MANE_114.groupby('Transcript stable ID version').agg({'Gene stable ID': 'first',
                                                                                         'cDNA coding start': 'min',
                                                                                         'cDNA coding end': 'max'})

Ensembl_canonical_starts['cDNA coding start'] = Ensembl_canonical_starts['cDNA coding start'].apply(
    lambda x: int(x) - 1)
Ensembl_canonical_starts['cDNA coding end'] = Ensembl_canonical_starts['cDNA coding end'].apply(
    lambda x: int(x))
Ensembl_canonical_starts = Ensembl_canonical_starts.reset_index()

Ensembl_canonical_starts.to_csv(f'{out_file_name}.csv')
Ensembl_canonical_starts.loc[:, ['Transcript stable ID version', 'cDNA coding start',
                                 'cDNA coding end']].to_csv(f'{out_file_name}.bed', sep='\t', index=False, header=False)
