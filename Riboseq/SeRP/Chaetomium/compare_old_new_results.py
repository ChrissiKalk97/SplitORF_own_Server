import sys
import pandas as pd


# chae_gtf_path = sys.argv[1]
# transcript_fasta = sys.argv[2]
# out_file_name = sys.argv[3]

old_results = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/Analysis_31_10_25_Remus/final_results_SND3_WT_LFC1_padj0.01.csv'
old_results_all = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/Analysis_31_10_25_Remus/final_results_SND3_WT_LFC1_padj0.01_all.csv'
new_results = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0_not_IP_WT_vs_IN_0.5_andpadj_0.05.txt'


E_S_over_E_WT = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/E_S_over_E_WT_1_0_all_results_p_smaller_0.01.csv'
E_over_In_S = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/E_over_In_S_1_0_all_results_p_smaller_0.01.csv'
E_over_In_WT = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/E_over_In_WT_1_0_all_results_p_smaller_0.01.csv'

original_suppl3_file = '/projects/serp/work/references/Supplementary_File-3.csv'

new_results_df = pd.read_csv(new_results, header=None)

old_results_df = pd.read_csv(old_results, header=0, sep=";")

old_results_all_df = pd.read_csv(old_results_all, header=0, sep=";")

# read in all diff regulated mRNAs

E_S_over_E_WT_df = pd.read_csv(E_S_over_E_WT, header=0)

E_over_In_S_df = pd.read_csv(E_over_In_S, header=0)

E_over_In_WT_df = pd.read_csv(E_over_In_WT, header=0)


# read in original suppl3 file
original_suppl3_file_df = pd.read_csv(original_suppl3_file, header=0, sep=";")

original_suppl3_file_filtered_df = original_suppl3_file_df[original_suppl3_file_df['Transcript ID'].isin(
    new_results_df.iloc[:, 0])]


E_over_In_S_log_dict = E_over_In_S_df.set_index(
    'Unnamed: 0')['log2FoldChange'].to_dict()
E_over_In_S_p_dict = E_over_In_S_df.set_index('Unnamed: 0')['padj'].to_dict()
original_suppl3_file_filtered_df['log2FoldChange_SND3vsInput'] = original_suppl3_file_filtered_df['Transcript ID'].map(
    E_over_In_S_log_dict)
original_suppl3_file_filtered_df['padj_SND3vsInput'] = original_suppl3_file_filtered_df['Transcript ID'].map(
    E_over_In_S_p_dict)


E_S_over_E_WT_log_dict = E_S_over_E_WT_df.set_index(
    'Unnamed: 0')['log2FoldChange'].to_dict()
E_S_over_E_WT_p_dict = E_S_over_E_WT_df.set_index('Unnamed: 0')[
    'padj'].to_dict()
original_suppl3_file_filtered_df['log2FoldChange_SND3vsWT'] = original_suppl3_file_filtered_df['Transcript ID'].map(
    E_S_over_E_WT_log_dict)
original_suppl3_file_filtered_df['padj_SND3vsWT'] = original_suppl3_file_filtered_df['Transcript ID'].map(
    E_S_over_E_WT_p_dict)


E_over_In_WT_log_dict = E_over_In_WT_df.set_index(
    'Unnamed: 0')['log2FoldChange'].to_dict()
E_over_In_WT_p_dict = E_over_In_WT_df.set_index('Unnamed: 0')['padj'].to_dict()
original_suppl3_file_filtered_df['log2FoldChange_WTvsInput'] = original_suppl3_file_filtered_df['Transcript ID'].map(
    E_over_In_WT_log_dict)
original_suppl3_file_filtered_df['padj_WTvsInput'] = original_suppl3_file_filtered_df['Transcript ID'].map(
    E_over_In_WT_p_dict)


print('Number of transcripts that are found in old and new analysis:', sum(
    old_results_df['Transcript ID'].isin(original_suppl3_file_filtered_df['Transcript ID'])))

assert (
    original_suppl3_file_filtered_df['log2FoldChange_SND3vsInput'] > 1.0).all()
assert (
    original_suppl3_file_filtered_df['log2FoldChange_SND3vsWT'] > 1.0).all()
assert (
    (original_suppl3_file_filtered_df['log2FoldChange_WTvsInput'] < 1.0) | (original_suppl3_file_filtered_df['padj_WTvsInput'] > 0.05) | (original_suppl3_file_filtered_df['log2FoldChange_WTvsInput'].isna())).all()

assert (original_suppl3_file_filtered_df['padj_SND3vsInput'] < 0.01).all()
assert (original_suppl3_file_filtered_df['padj_SND3vsWT'] < 0.01).all()

original_suppl3_file_filtered_df.to_csv(
    '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0_not_IP_WT_vs_IN_0.5_andpadj_0.05_with_information.csv')
