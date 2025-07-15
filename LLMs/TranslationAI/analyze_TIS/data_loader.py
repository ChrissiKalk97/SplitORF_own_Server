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


def load_so_results(so_results):
    predicted_so_orfs = pd.read_csv(so_results, header=0, sep='\t')
    so_transcripts = predicted_so_orfs['OrfTransID'].to_list()
    predicted_so_orfs['OrfPos'] = predicted_so_orfs['OrfPos'].apply(
        lambda x: x.split(','))
    predicted_so_orfs['OrfStarts'] = predicted_so_orfs['OrfPos'].apply(
        lambda x: [y.split('-')[0] for y in x])
    # get Id to match the TIS transformer output one!
    predicted_so_orfs['nr_SO_starts'] = predicted_so_orfs['OrfPos'].apply(
        lambda x: len(x))
    return predicted_so_orfs, so_transcripts


def load_Ensembl_canonical(Ensembl_path):
    Ensembl_canonical_starts = pd.read_csv(
        Ensembl_path, sep='\t', header=0, encoding="utf-8")
    Ensembl_canonical_starts = Ensembl_canonical_starts.groupby('Transcript stable ID').agg({'Gene stable ID': 'first',
                                                                                             'cDNA coding start': 'min',
                                                                                             'cDNA coding end': 'max'})
    Ensembl_canonical_starts['cDNA coding start'] = Ensembl_canonical_starts['cDNA coding start'].apply(
        lambda x: int(x) - 1)
    Ensembl_canonical_starts['cDNA coding end'] = Ensembl_canonical_starts['cDNA coding end'].apply(
        lambda x: int(x))
    Ensembl_canonical_starts = Ensembl_canonical_starts.reset_index()
    Ensembl_canonical_starts.rename(
        columns={'Transcript stable ID': 'OrfTransID'}, inplace=True)
    return Ensembl_canonical_starts


def load_dna_ur_df(UR_path):
    dna_ur_df = pd.read_csv(UR_path, sep='\t', header=None, names=[
                            'chr', 'start', 'stop', 'ID', 'score', 'strand'])
    dna_ur_df['OrfID'] = dna_ur_df['ID'].str.split(':').apply(lambda x: x[1])
    dna_ur_df['OrfTransID'] = dna_ur_df['ID'].str.split(
        ':').apply(lambda x: x[0])

    nr_orfs_with_UR = len(dna_ur_df['OrfID'].unique())
    nr_transcripts_with_UR = len(dna_ur_df['OrfTransID'].unique())
    print('Number of ORFs with unique region', nr_orfs_with_UR)
    print('Number of transcripts with unique region', nr_transcripts_with_UR)

    return dna_ur_df, nr_orfs_with_UR, nr_transcripts_with_UR
