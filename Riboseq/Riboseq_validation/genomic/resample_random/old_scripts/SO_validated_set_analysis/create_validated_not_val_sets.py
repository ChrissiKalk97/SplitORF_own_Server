import sys
import os
import re
import glob

import seaborn as sbn
import matplotlib.pyplot as plt
import pandas as pd


sys.path.append('/home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/analyze_TIS')

from data_loader import load_TranslationAI, load_SO_results  # noqa: F401
from analysis import get_TransAI_SO_info, get_orf_start_probs  # noqa: F401
from plotting import plot_emp_background_TransAI  # noqa: F401

SO_results = '/projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt'
Ribo_coverage_path = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/NMD_genome'
DNA_UR_path = '/projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_NMD_16_12_24.bed'
out_dir = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/SO_validated_set_analysis'
region_type = 'NMD'
TIS_results = '/projects/splitorfs/work/LLMs/TranslationAI/Output/NMD_trnascripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt'


################################################################################
##
# Helper Functions
##
################################################################################

def plot_three_category_pie(cat1, cat2, cat3, nr_total_SO, names_list, title, out_dir, figname, region_type, color_order):
    # plot pie chart with numbers: first, middle, last ORF
    fig, ax = plt.subplots()
    ax.pie([cat1, cat2, cat3],
           labels=names_list,
           autopct=lambda p: '{:.0f}'.format(p * nr_total_SO / 100),
           colors=color_order)
    plt.title(f'{title}')

    plt.savefig(
        os.path.join(out_dir, f'{region_type}_{figname}.png'))
    plt.close()


def get_SO_position_in_transcript(SO_df):
    # sort the ORF starts by position
    SO_df['OrfStarts'] = SO_df.apply(
        lambda x: sorted([int(start) for start in x['OrfStarts']]), axis=1)
    # map the ORF start to the respective position in the sorted list
    # indicate whether it is the first or a later (first, middle, last)
    SO_df['OrfIndex'] = SO_df.apply(
        lambda x: x['OrfStarts'].index(int(x['OrfStart'])), axis=1)
    SO_df['OrfPosition'] = SO_df.apply(lambda x: 'first' if x['OrfIndex'] == 0 else (
        'last' if x['OrfIndex'] == len(x['OrfStarts'])-1 else 'middle'), axis=1)
    return SO_df


def count_orfs_by_position(SO_df):
    nr_first_ORFs = sum(SO_df['OrfPosition'] == 'first')
    nr_middle_ORFs = sum(SO_df['OrfPosition'] == 'middle')
    nr_last_ORFs = sum(SO_df['OrfPosition'] == 'last')
    return nr_first_ORFs, nr_middle_ORFs, nr_last_ORFs


def plot_start_prob_by_orf_position(so_df, title, region_type, out_dir, hue='OrfPosition', palette=['#75C1C5', '#CC79A7', '#FFC500']):
    plt.figure(figsize=(8, 6))
    sbn.histplot(data=so_df, x="OrfProb",
                 hue=hue, palette=palette, multiple='stack')

    # Labels and legend
    plt.xlabel("Probability of SO start sites")
    plt.ylabel("Frequency")
    plt.title(
        f"{title} - {region_type}")

    plt.savefig(f"{out_dir}/{region_type}_Riboseq_empirical_SO_start_probs_by_position.svg",
                dpi=300, bbox_inches="tight")
    plt.close()


################################################################################
##
# Load all predicted split-ORFs
##
################################################################################
predicted_SO_ORFs, SO_transcripts = load_SO_results(SO_results)
predicted_SO_ORFs = predicted_SO_ORFs[['OrfTransID', 'OrfIDs', 'OrfStarts']]
predicted_SO_ORFs['OrfID'] = predicted_SO_ORFs.apply(
    lambda x: x['OrfIDs'].split(','), axis=1)
predicted_SO_ORFs['OrfStart'] = predicted_SO_ORFs['OrfStarts']
all_predicted_SO_ORFs = predicted_SO_ORFs.explode(
    ['OrfID', 'OrfStart'], ignore_index=True)

assert len(all_predicted_SO_ORFs['OrfTransID'].unique()) == 9261 or len(
    all_predicted_SO_ORFs['OrfTransID'].unique()) == 5625


total_nr_SO = len(all_predicted_SO_ORFs.index)
total_nr_SO_transcripts = len(all_predicted_SO_ORFs['OrfTransID'].unique())
print('Total number of predicted Split-ORFs', total_nr_SO)
print('Total number of predicted Split-ORF transcripts',
      total_nr_SO_transcripts)

# ------------------ LOAD VALIDATED SOs ------------------ #
for empirical_Ribo_findings_file in glob.glob(f"{Ribo_coverage_path}/*_unique_regions.csv"):
    print(empirical_Ribo_findings_file)
    sample = os.path.basename(
        empirical_Ribo_findings_file).rsplit('_', 2)[0]
    print(sample)

    all_predicted_SO_ORFs[sample] = 0

    Ribo_results_df = pd.read_csv(
        empirical_Ribo_findings_file, header=0, index_col=0)

    Ribo_results_significant_df = Ribo_results_df[Ribo_results_df["significant"] == 1]
    Ribo_results_significant_df['ORF'] = Ribo_results_significant_df['name'].apply(
        lambda x: x.split(':')[1])
    Ribo_results_significant_df['OrfTransID'] = Ribo_results_significant_df['name'].apply(
        lambda x: x.split(':')[0].split('|')[1])

    all_predicted_SO_ORFs[sample] = all_predicted_SO_ORFs['OrfID'].isin(
        Ribo_results_significant_df['ORF'])

# ------------ FILTER VAL SOs FOR VALIDATION IN 2 OR MORE SAMPLES -------------#
all_predicted_SO_ORFs['ValidationCount'] = all_predicted_SO_ORFs.iloc[:, range(
    5, 23)].sum(axis=1, numeric_only=True)

validated_SO_df = all_predicted_SO_ORFs[all_predicted_SO_ORFs['ValidationCount'] > 1]

nr_validated_SO = len(validated_SO_df.index)
nr_validated_transcripts = len(validated_SO_df['OrfTransID'].unique())

print('Number of validated SOs:', nr_validated_SO)
print('Number of validated SO transcripts:', nr_validated_transcripts)

# ------------------ LOAD DNA UNIQUE REGIONS ------------------ #
DNA_UR_df = pd.read_csv(DNA_UR_path, sep='\t', header=None, names=[
                        'chr', 'start', 'stop', 'ID', 'score', 'strand'])
DNA_UR_df['OrfID'] = DNA_UR_df['ID'].str.split(':').apply(lambda x: x[1])
DNA_UR_df['OrfTransID'] = DNA_UR_df['ID'].str.split(':').apply(lambda x: x[0])

nr_ORFs_with_UR = len(DNA_UR_df['OrfID'].unique())
nr_transcripts_with_UR = len(DNA_UR_df['OrfTransID'].unique())
print('Number of ORFs with unique region', nr_ORFs_with_UR)
print('Number of transcripts with unique region', nr_transcripts_with_UR)

# ------------------ PLOT VALIDATED TRANSCRIPT PROPORTIONS ------------------ #
nr_SO_not_validated = nr_ORFs_with_UR - nr_validated_SO
nr_SO_no_UR = total_nr_SO - nr_ORFs_with_UR

assert nr_validated_SO + nr_SO_not_validated + nr_SO_no_UR == total_nr_SO


plot_three_category_pie(nr_validated_SO,
                        nr_SO_not_validated,
                        nr_SO_no_UR,
                        total_nr_SO,
                        ['# validated SO', '# SO with UR not validated',
                            '# SO without UR'],
                        'Split-ORF validation pie chart',
                        out_dir,
                        'SO_validation_pie_chart',
                        region_type,
                        ['#CC79A7', '#FFC500', '#75C1C5']
                        )

# ------------------ PLOT VALIDATED SO POSITIONS ------------------ #
validated_SO_df = get_SO_position_in_transcript(validated_SO_df)
nr_first_ORFs, nr_middle_ORFs, nr_last_ORFs = count_orfs_by_position(
    validated_SO_df)

assert nr_first_ORFs + nr_middle_ORFs + nr_last_ORFs == nr_validated_SO

plot_three_category_pie(nr_first_ORFs,
                        nr_middle_ORFs,
                        nr_last_ORFs,
                        nr_validated_SO,
                        ['# first ORFs', '# middle ORFs', '# last ORFs'],
                        'Pie chart of ORF positions of validated Split-ORFs',
                        out_dir,
                        'positions_of_validated_SO_pie_chart',
                        region_type,
                        ['#75C1C5', '#CC79A7', '#FFC500']
                        )


################################################################################
##
# Plot distribution of DNA URs among first, middle and last ORFs
##
################################################################################
DNA_UR_df = pd.merge(
    DNA_UR_df, all_predicted_SO_ORFs.iloc[:, 1:5], on='OrfID', how='left')

nr_DNA_URs = len(DNA_UR_df.index)

DNA_UR_df = get_SO_position_in_transcript(DNA_UR_df)
nr_first_ORFs, nr_middle_ORFs, nr_last_ORFs = count_orfs_by_position(
    DNA_UR_df)
plot_three_category_pie(nr_first_ORFs,
                        nr_middle_ORFs,
                        nr_last_ORFs,
                        nr_DNA_URs,
                        ['# first ORFs', '# middle ORFs', '# last ORFs'],
                        'Pie chart of ORF positions of SOs with DNA UR',
                        out_dir,
                        'positions_of_SO_with_UR_pie_chart',
                        region_type,
                        ['#75C1C5', '#CC79A7', '#FFC500']
                        )


# ------------------ PLOT VALIDATED SO TRANSAI PROBS ------------------ #
trans_ai_df = load_TranslationAI(TIS_results)
trans_ai_df = get_TransAI_SO_info(trans_ai_df)
trans_ai_df = trans_ai_df.loc[:, ['OrfTransID', 'TIS_dict', 'TIS_pos_list']]
validated_SO_df_subset = validated_SO_df.loc[:, [
    'OrfTransID', 'OrfIDs', 'OrfStarts', 'OrfStart', 'OrfIndex', 'OrfPosition']]
validated_SO_df_subset['OrfNumber'] = validated_SO_df_subset['OrfIndex'] + 1


val_so_trans_ai_df, Ribo_df_merged_first, Ribo_df_merged_second = get_orf_start_probs(
    validated_SO_df_subset, trans_ai_df, on='OrfTransID')


plot_start_prob_by_orf_position(val_so_trans_ai_df,
                                'Probability of Riboseq validated ORFs by ORF position',
                                region_type,
                                out_dir)

plot_emp_background_TransAI(
    Ribo_df_merged_first, Ribo_df_merged_second, out_dir, region_type, 'all_validated')
