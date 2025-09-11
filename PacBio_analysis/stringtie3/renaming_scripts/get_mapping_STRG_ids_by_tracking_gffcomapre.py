#This script takes a gtf of an assembly from stringtie that was merged of single sample assemblies
#and writes the ref_gene_id attribute as the gene_id
import sys
import regex as re
import pandas as pd


colnames = ['new_tid', 'new_gid', 'original_name', 'similarity', 'gtf1', 'gtf2']
tracking_df = pd.read_csv(sys.argv[1], sep = '\t', names = colnames)
#print(tracking_df.head())

##########################################################################################
# FIRST GTF FILE                                                                        ##
##########################################################################################

#get gene id mapping
filtered_df = tracking_df[tracking_df.loc[:, 'gtf1'].str.split('|').str[0].str.contains('STRG')].copy()
filtered_df['STRG_gene_name'] = filtered_df['gtf1'].str.split('|').str[0].str.split(':').str[1]
filtered_df.loc[:, 'STRG_gene_name']  = filtered_df.loc[:,'gtf1'].str.split('|').str[0].str.split(':').str[1]
gene_dict_gtf1 = dict(zip(filtered_df['STRG_gene_name'], filtered_df['new_gid']))
#print(gene_dict_gtf1['STRG.19499'])

#get transcript id mapping
filtered_tid_df = tracking_df[tracking_df.loc[:, 'gtf1'].str.split('|').apply(len) > 1].copy()
filtered_tid_df = filtered_tid_df[filtered_tid_df.loc[:, 'gtf1'].str.split('|').str[1].str.contains('STRG')]
filtered_tid_df['STRG_tid'] = filtered_tid_df.loc[:, 'gtf1'].str.split('|').str[1]
tid_dict_gtf1 = dict(zip(filtered_tid_df['STRG_tid'], filtered_tid_df['new_tid']))
#print(tid_dict_gtf1)



##########################################################################################
# SECOND GTF FILE                                                                        ##
##########################################################################################

#get gene id mapping
filtered_df2 = tracking_df[tracking_df.loc[:, 'gtf2'].str.split('|').str[0].str.contains('STRG')].copy()
filtered_df2.loc[:, 'STRG_gene_name'] = filtered_df2.loc[:, 'gtf2'].str.split('|').str[0].str.split(':').str[1]
gene_dict_gtf2 = dict(zip(filtered_df2['STRG_gene_name'], filtered_df2['new_gid']))
#print(gene_dict_gtf1)

#get transcript id mapping
filtered_tid_df2 = tracking_df[tracking_df.loc[:, 'gtf2'].str.split('|').apply(len) > 1].copy()
filtered_tid_df2 = filtered_tid_df2[filtered_tid_df2.loc[:, 'gtf2'].str.split('|').str[1].str.contains('STRG')]
filtered_tid_df2['STRG_tid'] = filtered_tid_df2.loc[:, 'gtf2'].str.split('|').str[1]
tid_dict_gtf2 = dict(zip(filtered_tid_df2['STRG_tid'], filtered_tid_df2['new_tid']))
#print(tid_dict_gtf1)



################################################################################
# WRITE FIRST MAPPED GTF FILE                                                 ##
################################################################################
assembled_gtf1 = open(sys.argv[2], 'r')
Lines = assembled_gtf1.readlines()
gtf_with_repalcement1 = open(sys.argv[3], 'w')

replaced_lines = []


for line in Lines:
    ref_gene_id = re.search(r'gene_id "(STRG\.\d*)"', line)
    ref_tid = re.search(r'transcript_id "(STRG\.\d*\.\d+)"', line)
    if ref_gene_id is not None:
        try:
            line = re.sub(r'"STRG\.\d*"', '"'+gene_dict_gtf1[ref_gene_id.groups(1)[0]]+'"', line)
        except:
            print(ref_gene_id.groups(1)[0])
        if ref_tid is not None:
            try:
                line = re.sub(r'"STRG\.\d+\.\d+"', '"'+tid_dict_gtf1[ref_tid.groups(1)[0]]+'"', line)
            except:
                print(ref_tid.groups(1)[0])
        replaced_lines.append(line)
    else:
        replaced_lines.append(line)

gtf_with_repalcement1.writelines(replaced_lines)



##########################################################################################
# WRITE SECOND MAPPED GTF FILE                                                                        ##
##########################################################################################
assembled_gtf2 = open(sys.argv[4], 'r')
Lines2 = assembled_gtf2.readlines()
gtf_with_repalcement2 = open(sys.argv[5], 'w')

replaced_lines2 = []


for line in Lines2:
    ref_gene_id = re.search(r'gene_id "(STRG\.\d*)"', line)
    ref_tid = re.search(r'transcript_id "(STRG\.\d*\.\d+)"', line)
    if ref_gene_id is not None:
        try:
            line = re.sub(r'"STRG\.\d*"', '"'+gene_dict_gtf2[ref_gene_id.groups(1)[0]]+'"', line)
        except:
            print('gtf2', ref_gene_id.groups(1)[0])
        if ref_tid is not None:
            try:
                line = re.sub(r'"STRG\.\d+\.\d+"', '"'+tid_dict_gtf2[ref_tid.groups(1)[0]]+'"', line)
            except:
                print('gtf2', ref_tid.groups(1)[0])
        replaced_lines2.append(line)
    else:
        replaced_lines2.append(line)

gtf_with_repalcement2.writelines(replaced_lines2)