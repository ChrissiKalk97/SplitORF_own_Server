import sys
import pandas as pd

UR_path = sys.argv[1]
outname = sys.argv[2]
# UR_path = "/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_16.12.2024-14.13.38_NMD_genomic_correction/Unique_DNA_regions_genomic.bed"

UR_df = pd.read_csv(UR_path, sep='\t', header=None, names=[
                    'chr', 'start', 'stop', 'ORF_ID', 'score', 'strand'])

# So I think that the start is independent of the strand and I always need
# to increase the start by one to get from bed to Ensembl format
UR_df['coordinate'] = (UR_df['start'] + 1).astype(
    str) + '-' + UR_df['stop'].astype(str)

UR_grouped_df = UR_df.groupby('ORF_ID').agg({
    'coordinate': ','.join})
UR_grouped_df['chrom'] = UR_df.groupby('ORF_ID').agg({'chr': 'first'})
UR_grouped_df['strand'] = UR_df.groupby('ORF_ID').agg({'strand': 'first'})
UR_grouped_df['ORF_type'] = 'Split-ORF'
UR_grouped_df['transcript_type'] = 'Split-ORF_transcript'

UR_grouped_df = UR_grouped_df.reset_index()
UR_grouped_df['transcript_id'] = UR_grouped_df['ORF_ID'].apply(
    lambda x: x.split('|')[1].split(':')[0])
UR_grouped_df['gene_id'] = UR_grouped_df['ORF_ID'].apply(
    lambda x: x.split('|')[0])
UR_grouped_df['gene_name'] = 'blablabla'
UR_grouped_df['gene_type'] = 'blublublu'
UR_grouped_df['start_codon'] = 'NNN'

UR_grouped_df = UR_grouped_df.loc[:, ['ORF_ID', 'ORF_type', 'transcript_id', 'transcript_type',
                                      'gene_id', 'gene_name', 'gene_type', 'chrom', 'strand', 'start_codon', 'coordinate']]

UR_grouped_df.to_csv(outname, sep='\t', index=False)
