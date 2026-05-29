
import sys
from pygtftk.gtf_interface import GTF
from typing import List

def get_gene_string(gene_ids: List[str]) -> str:
    '''get string of transcript or other ids for filtering 
    of a GTF class object from pygtftk'''
    gene_string = ''
    if len(gene_ids) > 1:
        for gene_id in gene_ids[:-1]:
            gene_string += gene_id+','
    gene_string += gene_ids[-1]
    return gene_string


reference_gtf = GTF(sys.argv[1], check_ensembl_format=False)
custom_gtf = GTF(sys.argv[2], check_ensembl_format=False)

gene_ids = custom_gtf.get_gn_ids(nr=True)
gene_string = get_gene_string(gene_ids)

reference_gtf_CDS = reference_gtf\
.select_by_key('feature', 'CDS,exon,start_codon,stop_codon,transcript,gene')\
.select_by_key('transcript_biotype', 'protein_coding')\
.select_by_key('gene_id', gene_string)


reference_gtf_CDS.write(sys.argv[3])