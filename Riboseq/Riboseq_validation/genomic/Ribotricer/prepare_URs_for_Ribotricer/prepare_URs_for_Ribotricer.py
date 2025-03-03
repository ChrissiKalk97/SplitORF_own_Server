import sys
import pandas as pd

sys.path.append("/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/RiboTISH") 
from prepare_unique_positions_as_ORF_input import create_region_for_RiboTISH

UR_path = sys.argv[1]
outname = sys.argv[2]
# UR_path = "/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_16.12.2024-14.13.38_NMD_genomic_correction/Unique_DNA_regions_genomic.bed"

UR_df = pd.read_csv(UR_path, sep='\t', header=None, names=[
                    'chr', 'start', 'stop', 'ORF_ID', 'score', 'strand'])


UR_df['trans_coords'] = UR_df['ORF_ID'].apply(lambda x: x.split(':')[2:])
UR_df['trans_coords_modified'] = UR_df['trans_coords'].apply(lambda x: create_region_for_RiboTISH(int(x[0]), int(x[2]), int(x[3])))

# mod start is smaller or equal to start, calculate the offset as start minus mod start
UR_df['start_offset'] = UR_df.apply(lambda x: int(x['trans_coords'][2]) - x['trans_coords_modified'][0], axis = 1)
# modified end is equal or greater to end, calculate offset as mod end minus end
UR_df['end_offset'] = UR_df.apply(lambda x: x['trans_coords_modified'][1] - int(x['trans_coords'][3]), axis = 1)


# So I think that the start is independent of the strand and I always need
# to increase the start by one to get from bed to Ensembl format
UR_df['coordinate'] = (UR_df['start'] + 1).astype(
    str) + '-' + UR_df['stop'].astype(str)
UR_df['length'] = UR_df['stop'] - UR_df['start']

UR_grouped_df = UR_df.groupby('ORF_ID').agg({
    'coordinate': ','.join})
UR_grouped_df['chrom'] = UR_df.groupby('ORF_ID').agg({'chr': 'first'})
UR_grouped_df['strand'] = UR_df.groupby('ORF_ID').agg({'strand': 'first'})
UR_grouped_df['ORF_type'] = 'novel'
UR_grouped_df['transcript_type'] = 'novel'

# keep the offset information
UR_grouped_df['start_offset'] = UR_df.groupby('ORF_ID').agg({'start_offset': 'first'})
UR_grouped_df['end_offset'] = UR_df.groupby('ORF_ID').agg({'end_offset': 'first'})

# get start stop and total length
UR_grouped_df['min_pos'] = UR_df.groupby('ORF_ID').agg({'start': 'min'})
UR_grouped_df['max_pos'] = UR_df.groupby('ORF_ID').agg({'stop': 'max'})
UR_grouped_df['length'] = UR_df.groupby('ORF_ID').agg({'length': 'sum'})


# adjust the start, stop and length to multiplicity of 3 and correct frame
def adjust_end_and_start(start_offset, min_position, end_offset, max_position, strand, coordinates, length):
    if strand == '+':
        max_position = max_position + end_offset
        # min posiition is bed format, they require Ensembl
        min_position = min_position - start_offset + 1
    else:
        max_position = max_position + start_offset
        # min posiition is bed format, they require Ensembl
        min_position = min_position - end_offset + 1
    if len(coordinates.split('-')) == 2:
        coordinates = '-'.join([str(min_position), str(max_position)])
    else:
        coordinates = coordinates.split(',')
        coordinates[0] = str(min_position) + '-' + coordinates[0].split('-')[1]
        coordinates[-1] = coordinates[-1].split('-')[0] + '-' + str(max_position)
        coordinates = ','.join(coordinates)
    
    length = length + end_offset + start_offset
    
    return min_position, max_position, coordinates, length

UR_grouped_df[['min_pos', 'max_pos', 'coordinate', 'length']] = UR_grouped_df.apply(
    lambda x: pd.Series(adjust_end_and_start(
        x['start_offset'], 
        x['min_pos'], 
        x['end_offset'], 
        x['max_pos'], 
        x['strand'], 
        x['coordinate'],
        x['length']
    )), 
    axis=1
)



UR_grouped_df = UR_grouped_df.reset_index()
UR_grouped_df['ORF_ID'] = UR_grouped_df['ORF_ID'] + '_' + (UR_grouped_df['min_pos'] + 1).astype(str) + '_' + UR_grouped_df['max_pos'].astype(str) + '_' + UR_grouped_df['length'].astype(str)

assert (UR_grouped_df['length'] % 3 == 0).all()

UR_grouped_df['transcript_id'] = UR_grouped_df['ORF_ID'].apply(
    lambda x: x.split('|')[1].split(':')[0])
UR_grouped_df['gene_id'] = UR_grouped_df['ORF_ID'].apply(
    lambda x: x.split('|')[0])
UR_grouped_df['gene_name'] = UR_grouped_df['gene_id']
UR_grouped_df['gene_type'] = 'novel'
UR_grouped_df['start_codon'] = 'NNN'

UR_grouped_df = UR_grouped_df.loc[:, ['ORF_ID', 'ORF_type', 'transcript_id', 'transcript_type',
                                      'gene_id', 'gene_name', 'gene_type', 'chrom', 'strand', 'start_codon', 'coordinate']]

UR_grouped_df.to_csv(outname, sep='\t', index=False)
