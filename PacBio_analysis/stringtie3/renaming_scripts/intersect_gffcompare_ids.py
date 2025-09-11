import sys
import regex as re
from pybedtools import BedTool  
import pandas as pd 
import pygtftk
from pygtftk.gtf_interface import GTF
def main():
    def change_MSTRG_gene_name(input_gtf_file, output_gtf_file, unique_tids, summed_counts):
        """helper function to write the mapped gene names into a new gtf file"""
        input_gtf = open(input_gtf_file,'r')
        output_gtf = open(output_gtf_file, 'w')
        for l in input_gtf: 
            if not l.startswith('#'):
                info = re.split(r'\s+|;', l)
                if info[12].strip('"') in unique_tids:
                    # Modify the gene name

                    l = re.sub(r'gene_id "(XLOC\_[0-9]*)"','gene_id "'+ summed_counts[summed_counts['tid'] == info[12].strip('"')]['name'].values[0] + '"',l)
                # Write the feature to the output GTF file
            output_gtf.write(l)



    reference_bed = BedTool(sys.argv[2])
    MSTRG_transcripts_bed = BedTool(sys.argv[1])
    intersection = MSTRG_transcripts_bed.intersect(reference_bed, wao = True, s = True).saveas('intersection.bed')
    summed_counts = pd.read_table(intersection.fn, names=['chrom', 'start', 'stop', 'name', 'score', 'strand',\
                    'chrom_tar', 'start_tar', 'stop_tar', 'name_tar', 'score_tar', 'strand_tar', 'overlap'])
    print(summed_counts.head())
    summed_counts = summed_counts.groupby(['name', 'name_tar'])['overlap'].sum()
    summed_counts = summed_counts.reset_index()
    summed_counts =  summed_counts[summed_counts['overlap'] > 0]
    print(summed_counts.head())
    
    summed_counts['tid'] = summed_counts['name_tar'].str.split(':').str[1]
    rowIds = summed_counts.groupby('tid')['overlap'].idxmax()
    summed_counts = summed_counts.loc[rowIds]
    unique_tids = summed_counts['tid'].unique()

    input_gtf_file = sys.argv[3]
    output_gtf_file = sys.argv[4]
    
    change_MSTRG_gene_name(input_gtf_file, output_gtf_file, unique_tids, summed_counts)



if __name__ == "__main__":
    main()