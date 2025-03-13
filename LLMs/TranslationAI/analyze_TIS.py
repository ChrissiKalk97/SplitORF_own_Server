import sys
import pandas as pd
import re

TIS_results = sys.argv[1]
SO_results = sys.argv[2]

# TIS_results = '/projects/splitorfs/work/LLMs/TranslationAI/Output/NMD_trnascripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt'
# SO_results = '/projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt'

# only read the first 10 start codons, the rest will not be of interest
TIS_results_df = pd.read_csv(TIS_results, sep="\t", usecols=range(11), header = None)
TIS_results_df.rename(columns={0:'fasta_header'}, inplace=True)
TIS_results_df['OrfTransID'] = TIS_results_df['fasta_header'].apply(lambda x: re.split(r'\(|\)|\|', x)[4])


predicted_SO_ORFs = pd.read_csv(SO_results, header=0, sep='\t')
SO_transcripts = predicted_SO_ORFs['OrfTransID'].to_list()
predicted_SO_ORFs['OrfPos'] = predicted_SO_ORFs['OrfPos'].apply(
    lambda x: x.split(','))
predicted_SO_ORFs['OrfStarts'] = predicted_SO_ORFs['OrfPos'].apply(
    lambda x: [y.split('-')[0] for y in x])
# get Id to match the TIS transformer output one!
predicted_SO_ORFs['nr_SO_starts'] = predicted_SO_ORFs['OrfPos'].apply(lambda x: len(x))


TransAI_SO_preds = TIS_results_df[TIS_results_df['OrfTransID'].isin(SO_transcripts)]
print('Number of SO predicted transcripts', len(TransAI_SO_preds.index))
print('Number of non-SO predicted transcripts', len(TIS_results_df.index) - len(TransAI_SO_preds.index))

SO_2nd_start_mean_prob = TIS_results_df[TIS_results_df['OrfTransID'].isin(SO_transcripts)].iloc[:,2].dropna().apply(lambda x: float(x.split(',')[1])).mean()

random_2nd_start_mean_prob = TIS_results_df[~TIS_results_df['OrfTransID'].isin(SO_transcripts)].iloc[:,2].dropna().apply(lambda x: float(x.split(',')[1])).mean()

print('Mean probability of 2nd start in SO transcripts', SO_2nd_start_mean_prob)
print('Mean probability of 2nd starts in non-SO predicted NMD transcripts', random_2nd_start_mean_prob)


# How often do positions corresponds?
TransAI_SO_preds["TIS_dict"] = TransAI_SO_preds[range(1,11)].apply(lambda row: {str(k).split(',')[0]: str(k).split(',')[1] for k in row.dropna()}, axis=1)
TransAI_SO_preds["TIS_pos_list"] = TransAI_SO_preds[range(1,11)].apply(lambda row: [str(k).split(',')[0] for k in row.dropna()], axis=1)


df_merged = pd.merge(TransAI_SO_preds[["TIS_dict", "TIS_pos_list", "OrfTransID"]], predicted_SO_ORFs, on='OrfTransID', how='outer')

# Are k ORFs predicted in the top k?
# First sample the TranslationAI predicted positions for the same number as the predicted SOs for that transcript
df_merged['top_k_TIS'] = df_merged.apply(lambda x: x['TIS_pos_list'][:len(x['OrfStarts'])], axis = 1)
df_merged['SOs_are_top_TIS'] = df_merged.apply(lambda x: all(start in x['top_k_TIS'] for start in x['OrfStarts']), axis = 1)

print('Number of times that all k SO starts are within the topk predicted TIS by TranslationaAI', df_merged['SOs_are_top_TIS'].sum())

# Find out how often at least 2 predicted SO start positions obtain a probability of higher than or equal to 10%
# this will get the predicted probabilities for the SO start sites in the order that they are listed in 'OrfStarts'
df_merged['OrfStartsTransAIProbs'] = df_merged.apply(lambda x: [x['TIS_dict'][start] if start in x['TIS_pos_list'] else 0 for start in x['OrfStarts']], axis = 1)
df_merged['NrOrfStartsGreater0_1'] = df_merged['OrfStartsTransAIProbs'].apply(lambda x: len([prob for prob in x if float(prob) >= 0.1]))
print("Number of times at least 2 start codons of the SO pipeline have at least 0.1 prediction probability by TranslationAI: ", sum(df_merged['NrOrfStartsGreater0_1'] > 1))

