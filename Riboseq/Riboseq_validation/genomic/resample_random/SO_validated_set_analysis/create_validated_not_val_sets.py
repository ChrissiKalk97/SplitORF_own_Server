
import sys
import os
import re
import glob

import seaborn as sbn
import matplotlib.pyplot as plt
import pandas as pd


sys.path.append('/home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/analyze_TIS')

from data_loader import load_TranslationAI, load_SO_results  # noqa: F401


SO_results = '/projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt'
Ribo_coverage_path = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/NMD_genome'
DNA_UR_path = '/projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_NMD_16_12_24.bed'
out_dir = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/SO_validated_set_analysis'
region_type = 'NMD'

# ------------------ LOAD ALL PREDICTED SOs ------------------ #
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
DNA_UR_df['ORF_nr'] = DNA_UR_df['ID'].str.split(':').apply(lambda x: x[1])
DNA_UR_df['OrfTransID'] = DNA_UR_df['ID'].str.split(':').apply(lambda x: x[0])

nr_ORFs_with_UR = len(DNA_UR_df['ORF_nr'].unique())
nr_transcripts_with_UR = len(DNA_UR_df['OrfTransID'].unique())
print('Number of ORFs with unique region', nr_ORFs_with_UR)
print('Number of transcripts with unique region', nr_transcripts_with_UR)

# ------------------ PLOT VALIDATED TRANSCRIPT PROPORTIONS ------------------ #
nr_SO_not_validated = nr_ORFs_with_UR - nr_validated_SO
nr_SO_no_UR = total_nr_SO - nr_ORFs_with_UR

fig, ax = plt.subplots()
ax.pie([nr_validated_SO, nr_SO_not_validated, nr_SO_no_UR],
       labels=['# validated SO', '# SO with UR not validated', '# SO without UR'],
       autopct=lambda p: '{:.0f}'.format(p * total_nr_SO / 100))
plt.title(f'Split-ORF validation pie chart')

plt.savefig(
    os.path.join(out_dir, f'{region_type}_SO_validation_pie_chart.png'))
plt.close()

# ------------------ PLOT VALIDATED SO POSITIONS ------------------ #
# sort the ORF starts by position
validated_SO_df['OrfStarts'] = validated_SO_df.apply(
    lambda x: sorted([int(start) for start in x['OrfStarts']]), axis=1)
# map the ORF start to the respective position in the sorted list
# indicate whether it is the first or a later (first, middle, last)
validated_SO_df['OrfIndex'] = validated_SO_df.apply(
    lambda x: x['OrfStarts'].index(int(x['OrfStart'])), axis=1)
validated_SO_df['OrfPosition'] = validated_SO_df.apply(lambda x: 'first' if x['OrfIndex'] == 0 else (
    'last' if x['OrfIndex'] == len(x['OrfStart'])-1 else 'middle'), axis=1)

nr_first_ORFs = sum(validated_SO_df['OrfPosition'] == 'first')
nr_middle_ORFs = sum(validated_SO_df['OrfPosition'] == 'middle')
nr_last_ORFs = sum(validated_SO_df['OrfPosition'] == 'last')

# plot pie chart with numbers: first, middle, last ORF
fig, ax = plt.subplots()
ax.pie([nr_first_ORFs, nr_middle_ORFs, nr_last_ORFs],
       labels=['# first ORFs', '# middle ORFs', '# last ORFs'],
       autopct=lambda p: '{:.0f}'.format(p * nr_validated_SO / 100))
plt.title(f'Pie chart of ORF positions of validated Split-ORFs')

plt.savefig(
    os.path.join(out_dir, f'{region_type}_positions_of_validated_SO_pie_chart.png'))
plt.close()


# ------------------ PLOT VALIDATED SO TRANSAI PROBS ------------------ #
