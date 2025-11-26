import pandas as pd


riboseq_org_metadata_csv = "RiboSeqOrg_Metadata.csv"

riboseq_org_metadata_df = pd.read_csv(riboseq_org_metadata_csv)

riboseq_org_metadata_df = riboseq_org_metadata_df[riboseq_org_metadata_df['avgLength'] > 28]

riboseq_org_metadata_df = riboseq_org_metadata_df[riboseq_org_metadata_df['LIBRARYTYPE'].str.lower(
) == 'ribo-seq']


bioproject_selection_df = riboseq_org_metadata_df.groupby(
    'BioProject').agg({'Run': 'count'}) > 9


bioprojects_chosen = bioproject_selection_df[bioproject_selection_df['Run']].index

riboseq_org_metadata_df_filtered = riboseq_org_metadata_df[riboseq_org_metadata_df['BioProject'].isin(
    bioprojects_chosen.to_list())]

set(riboseq_org_metadata_df_filtered['BioProject'])

# there would be the possibilty to only consider data with UMIs
riboseq_org_metadata_df_filtered[riboseq_org_metadata_df_filtered['UMI']
                                 != '0.0']['BioProject']
len(set(
    riboseq_org_metadata_df_filtered[riboseq_org_metadata_df_filtered['UMI'] != '0.0']['BioProject']))
# 58 rows
# 2 projects, PRJNA725118 is mitochondrion profiling
# PRJNA928376, the others are actually the Ingolia data

# wha tif we do not require a certain number of samples
set(riboseq_org_metadata_df[riboseq_org_metadata_df['UMI']
    != '0.0']['BioProject'])
# still only these two...


# Maybe rather select by cell type
cell_type_list = list(set(riboseq_org_metadata_df_filtered['CELL_LINE']))
tissue_list = list(set(riboseq_org_metadata_df_filtered['TISSUE']))


# cancer Ribo-seq
# idea: could perform the analysis on all cancer datasets that are in the database
riboseq_org_metadata_df_filtered['Cancer'].unique()
# array(['0.0', 'no apparent disease', 'metastatic', 'primary tumor',
#        'fibrocystic disease', 'Glioblastoma,ud',
#        'activate B cell lymphoma (ABC)',
#        'germinal center B cell lymphoma (GCB)', 'poorly metastatic',
#        'highly metastatic to lungs', 'non-small cell lung carcinoma'],
#       dtype=object)
cancer_riboseq_df = riboseq_org_metadata_df_filtered[~riboseq_org_metadata_df_filtered['Cancer'].isin([
                                                                                                      '0.0', 'no apparent disease', 0.0])]

cancer_riboseq_df['BioProject'].unique()

# array(['PRJNA523167', 'PRJNA591767', 'PRJNA647736', 'PRJNA898352',
#        'PRJEB33244'], dtype=object)
