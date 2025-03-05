import sys
import os
import pandas as pd
import glob
import ast


FuFi_path = "/home/ckalk/scripts/FuFis/src"
sys.path.append(FuFi_path)

import GTF_Processing

gtf_path = "/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"


os.chdir("/projects/splitorfs/work/Riboseq/Output/Ribotricer/URs_as_ORFs")

def compare_reading_frames(coverage_list):
    RF1 = sum(SOs_filtered.loc[0,'profile'][::3])
    RF2 = sum(SOs_filtered.loc[0,'profile'][1::3])
    RF3 = sum(SOs_filtered.loc[0,'profile'][2::3])

    if RF1 > RF2 and RF1 > RF3:
        return True
    else:
        return False




for pred_ORFs_file in glob.glob("*_translating_ORFs.tsv"):

    sample_name = pred_ORFs_file.rsplit('_', maxsplit = 6)[0]
    print(sample_name)

    pred_ORFs_df = pd.read_csv(pred_ORFs_file, sep = "\t", header = 0)
    # select only the Sos
    SOs = pred_ORFs_df[pred_ORFs_df["ORF_type"] == "SO"]
    SOs['genomic_region'] = SOs['ORF_ID'].apply(lambda x: x.split('_', maxsplit = 1)[1])

    SOs_filtered = SOs[SOs["valid_codons_ratio"] > 0.3]
    SOs_filtered = SOs_filtered.reset_index()
    SOs_filtered['profile'] = SOs_filtered['profile'].apply(lambda x: ast.literal_eval(x))
    SOs_filtered['correct_reading_frame'] = SOs_filtered['profile'].apply(lambda x: compare_reading_frames(x))
    
    print(sample_name)
    print("Nr of predicted SOs with codon ratio higher 0.3: ", SOs_filtered['genomic_region'].nunique())
    print("Nr of predicted SOs no filter: ", SOs['genomic_region'].nunique())
    print('Nr URs with wrong frame', SOs_filtered[SOs_filtered['correct_reading_frame'] == False]['genomic_region'].nunique())
    # print(SOs_filtered.sort_values('phase_score'))


    # get the gene symbols out
    gene_symbol, _ = GTF_Processing.match_gene_identifiers(SOs_filtered['gene_id'].to_list(), gtf_file=gtf_path, species='human', scopes='ensembl.gene', fields='symbol', ensemblonly=False, return_full=False)
    gene_symbol_mapping = {k:v['symbol'] for (k,v) in gene_symbol.items()}
    SOs_filtered['gene_name'] = SOs_filtered['gene_id'].map(gene_symbol_mapping)

    # report only one entry per genomic region
    SOs_filtered = SOs_filtered.groupby('genomic_region').agg('first')
    SOs_filtered = SOs_filtered.reset_index()
    SOs_filtered[['ORF_ID', 'phase_score', 'read_count', 'valid_codons', 'valid_codons_ratio', 'gene_id', 'gene_name', 'chrom', 'correct_reading_frame']].to_csv(f'{sample_name}_URs_translating_0_3_cutoff_ribotricer.csv', index = False)
