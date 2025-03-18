import re
import pandas as pd


def load_TranslationAI(TIS_results):
    # only read the first 10 start codons, the rest will not be of interest
    TIS_results_df = pd.read_csv(
        TIS_results, sep="\t", usecols=range(11), header=None)
    TIS_results_df.rename(columns={0: 'fasta_header'}, inplace=True)
    TIS_results_df['OrfTransID'] = TIS_results_df['fasta_header'].apply(
        lambda x: re.split(r'\(|\)|\|', x)[4])

    return TIS_results_df


def load_SO_results(SO_results):
    predicted_SO_ORFs = pd.read_csv(SO_results, header=0, sep='\t')
    SO_transcripts = predicted_SO_ORFs['OrfTransID'].to_list()
    predicted_SO_ORFs['OrfPos'] = predicted_SO_ORFs['OrfPos'].apply(
        lambda x: x.split(','))
    predicted_SO_ORFs['OrfStarts'] = predicted_SO_ORFs['OrfPos'].apply(
        lambda x: [y.split('-')[0] for y in x])
    # get Id to match the TIS transformer output one!
    predicted_SO_ORFs['nr_SO_starts'] = predicted_SO_ORFs['OrfPos'].apply(
        lambda x: len(x))
    return predicted_SO_ORFs, SO_transcripts
