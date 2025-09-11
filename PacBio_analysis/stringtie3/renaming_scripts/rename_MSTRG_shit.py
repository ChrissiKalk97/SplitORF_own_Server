# This script takes a gtf of an assembly from stringtie that was merged of single sample assemblies
# and writes the ref_gene_id attribute as the gene_id
# and the ref_id as the transcript id
import sys
import regex as re
import pandas as pd
from pygtftk.gtf_interface import GTF

custom_gtf = open(sys.argv[1], 'r')
Lines = custom_gtf.readlines()
refseq_anno_file = pd.read_csv(sys.argv[2], sep="\t", comment='#', header = 0)
refseq_anno_file = refseq_anno_file[refseq_anno_file["ref_gene_id"] != "-"]
#print(refseq_anno_file[refseq_anno_file['ref_gene_id'].str.contains('-')])
gene_dict = pd.Series(refseq_anno_file.ref_gene_id.values,index=refseq_anno_file.qry_gene_id).to_dict()
#print(gene_dict)
refseq_equal_matches = refseq_anno_file[refseq_anno_file["class_code"] == "="]
refseq_equal_matches = refseq_equal_matches[refseq_equal_matches['qry_id'] != '-']
print(refseq_equal_matches.head())
tid_dict = pd.Series(refseq_equal_matches.ref_id.values,index=refseq_equal_matches.qry_id).to_dict()

gtf_with_repalcement = open(sys.argv[3], 'w')


replaced_lines = []

for line in Lines:
    if not line.startswith('#'):
        trans_id_match = re.search(r'transcript_id "(MSTRG\.\d*\.\d*)"', line)
        gene_id_match = re.search(r'gene_id "(MSTRG\.\d*)"', line)
        #print(gene_id_match.group(1)[0])
        if gene_id_match.group(1) in gene_dict.keys():
            if gene_dict[gene_id_match.group(1)] != '-':
                line = re.sub(r'gene_id "MSTRG\.\d*"', 'gene_id "'+gene_dict[gene_id_match.group(1)]+'"', line)
            if trans_id_match.group(1) in tid_dict.keys():
                if tid_dict[trans_id_match.group(1)] != '-':
                    line = re.sub(r'"MSTRG\.\d*\.\d*"', '"'+tid_dict[trans_id_match.group(1)]+'"', line)
    replaced_lines.append(line)
    

gtf_with_repalcement.writelines(replaced_lines)