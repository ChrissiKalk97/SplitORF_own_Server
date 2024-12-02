#get_CDS_coords.py is given structures information from Biomart on the CDS coordiantes
#required fields are: Gene ID, transcript ID, cDNA coding end, cDNA coding start and
#Transcript length (including UTRs and CDS)
#the exon wise information is combined such that the CDS start and end is extracted
#in transcript coordinates
#from this a bed file containing the 3' UTR regions is extracted for the random
#intersection in the Riboseq pipeline

import sys
import pandas as pd
from pybedtools import BedTool
from pygtftk.gtf_interface import GTF

#read in gtf file Ensembl that is tsl 1 and 2 filtered
Ensembl_tsl_filtered_gtf = GTF(sys.argv[1])
tsl_filtered_prot_cod_tids = Ensembl_tsl_filtered_gtf.select_by_key('transcript_support_level', '1').select_by_key('transcript_biotype', 'protein_coding').get_tx_ids(nr = True)


#read file with CDS cooridnates
Ensembl_cDNA_coding = pd.read_csv(sys.argv[2], sep='\t')
print(len(Ensembl_cDNA_coding['Transcript stable ID']))
#filter for tids protein coding and tsl 1and 2
Ensembl_cDNA_coding = Ensembl_cDNA_coding[Ensembl_cDNA_coding['Transcript stable ID'].isin(tsl_filtered_prot_cod_tids)]
print(len(Ensembl_cDNA_coding['Transcript stable ID']))
Ensembl_cDNA_coding = Ensembl_cDNA_coding[Ensembl_cDNA_coding['cDNA coding start'].notna()]
#create names for bed file
Ensembl_cDNA_coding['ID'] = Ensembl_cDNA_coding['Gene stable ID'] + '|' + Ensembl_cDNA_coding['Transcript stable ID'] 
Ensembl_cDNA_coding.drop(columns = ['Gene stable ID','Transcript stable ID'], inplace = True)
Ensembl_cDNA_coding.reset_index(inplace = True, drop = True)
#get the CDS start and end (min and max) of the exonwise information
start_idx = Ensembl_cDNA_coding.groupby('ID')['cDNA coding start'].idxmin().to_list()
stop_idx = Ensembl_cDNA_coding.groupby('ID')['cDNA coding end'].idxmax().to_list()
filteridx = start_idx + stop_idx

Ensembl_cDNA_coding = Ensembl_cDNA_coding.loc[filteridx]
Ensembl_cDNA_coding.reset_index(inplace = True)
Ensembl_cDNA_coding.sort_values(by = 'ID', axis = 0, inplace = True)
#print(Ensembl_cDNA_coding.head(50))

#combine entries to obtain CDS coordinates, min for start, max for end
d = {'cDNA coding end': 'max', 'cDNA coding start': 'min', 'CDS Length': 'first',\
      'Transcript length (including UTRs and CDS)': 'first'}
df_new = Ensembl_cDNA_coding.groupby('ID').aggregate(d).reset_index()

print(df_new.head)

#check is the CDS end = to start + length?
assert (df_new['cDNA coding end'] -  df_new['cDNA coding start'] + 1 == df_new['CDS Length']).all()

#df_new['3_prime_UTR'] = df_new['Transcript length (including UTRs and CDS)'] - df_new['cDNA coding end']
#df_new['5_prime_UTR'] = df_new['cDNA coding start']

#filter for transcript having either three prime of 5 prime
df_filtered = df_new[(df_new['cDNA coding end'] != df_new['Transcript length (including UTRs and CDS)'])]

#build bed string with 3' and 5' coordinates
coordinate_string = ''
for line in df_filtered.index:
    if df_filtered.loc[line, 'cDNA coding end'] < df_filtered.loc[line, 'Transcript length (including UTRs and CDS)']:
        coordinate_string += df_filtered.loc[line, 'ID']+'\t'+ str(int(df_filtered.loc[line, 'cDNA coding end']))+'\t'+\
        str(int(df_filtered.loc[line, 'Transcript length (including UTRs and CDS)']))+'\n'


coordinate_bed = BedTool(coordinate_string, from_string = True).saveas(sys.argv[3])
