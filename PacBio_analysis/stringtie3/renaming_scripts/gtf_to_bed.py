#gtf_to_bed.py takes a gtf provided as the first argument,
#filters it for the MSTRG entries  and writes it to a bed-6 file
#with the name provided by the second argument


import sys
from pybedtools import BedTool
from pygtftk.gtf_interface import GTF

def main():
    if sys.argv[3] == "MSTRG":
        gtf_file = BedTool(GTF(sys.argv[1], check_ensembl_format=False)\
            .select_by_key('feature', 'transcript')
            .to_bed(name=('gene_id', 'transcript_id'), sep=':'))
    
        gtf_file.filter(lambda gene:\
                        gene.name.split(":")[0].startswith('MSTRG')).saveas(sys.argv[2])
    elif sys.argv[3] == "reference":
        gtf_file = BedTool(GTF(sys.argv[1], check_ensembl_format=False)\
            .select_by_key('feature', 'gene')
            .to_bed(name='gene_id')).saveas(sys.argv[2])


if __name__ == "__main__":
    main()
