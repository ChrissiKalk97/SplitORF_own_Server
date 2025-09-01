import pandas as pd
classification_huvec = pd.read_csv(
    '/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/SQANTI3/SQANTI3_QC/HUVEC/isoforms_classification.txt', sep='\t')

classification_huvec_short = classification_huvec.loc[:, ['isoform',
                                                          'structural_category',
                                                          'diff_to_TSS',
                                                          'diff_to_TTS',
                                                          'diff_to_gene_TSS',
                                                          'diff_to_gene_TTS',
                                                          'subcategory',
                                                          'RTS_stage',
                                                          'all_canonical',
                                                          'min_sample_cov',
                                                          'min_cov',
                                                          'bite',
                                                          'iso_exp',
                                                          'gene_exp',
                                                          'ratio_exp',
                                                          'FSM_class',
                                                          'perc_A_downstream_TTS',
                                                          'ratio_TSS']]

classification_huvec[classification_huvec['bite']
                     == True]['structural_category'].value_counts()

classification_huvec[classification_huvec['RTS_stage']
                     == True]['structural_category'].value_counts()

classification_huvec[classification_huvec['min_cov']
                     < 2]['structural_category'].value_counts()

classification_huvec[(classification_huvec['perc_A_downstream_TTS']
                      > 60)]['structural_category'].value_counts()


classification_huvec[(classification_huvec['perc_A_downstream_TTS'] > 60) & (
    classification_huvec['ratio_TSS'] < 1.2)]['structural_category'].value_counts()

classification_huvec[(classification_huvec['perc_A_downstream_TTS'] > 60) | (
    classification_huvec['ratio_TSS'] < 1.2)]['structural_category'].value_counts()

classification_huvec[(classification_huvec['perc_A_downstream_TTS'] > 60) | (
    classification_huvec['ratio_TSS'] < 1.2)]['structural_category'].value_counts()


classification_huvec[(classification_huvec['perc_A_downstream_TTS'] > 60) | (
    classification_huvec['min_cov'] < 2)]['structural_category'].value_counts()


classification_huvec[(classification_huvec['perc_A_downstream_TTS'] > 60) | (
    classification_huvec['RTS_stage'] == True)]['structural_category'].value_counts()


classification_huvec.groupby('structural_category')['min_cov'].median()
classification_huvec.groupby('structural_category')[
    'perc_A_downstream_TTS'].median()


classification_huvec[classification_huvec['polyA_motif_found']
                     == False]['structural_category'].value_counts()

# filter 1: for the first rule applied to all: poly_A > 70 and no polyA motif
classification_huvec_filtered = classification_huvec[(
    classification_huvec['polyA_motif_found'] == True) | (classification_huvec['perc_A_downstream_TTS'] < 70)]


filter2 = (
    (classification_huvec_filtered['min_cov'] > 2) & (
        ((classification_huvec_filtered['RTS_stage'] == False) & (classification_huvec_filtered['bite'] == False)) |
        ((classification_huvec_filtered['RTS_stage'] == False) & (classification_huvec_filtered['ratio_TSS'] > 1.2)) |
        ((classification_huvec_filtered['bite'] == False) & (
            classification_huvec_filtered['ratio_TSS'] > 1.2))
    )
)
classification_huvec_filtered[
    filter2]['structural_category'].value_counts()


filter2_1 = (
    (classification_huvec_filtered['min_cov'] > 2) & (
        ((classification_huvec_filtered['RTS_stage'] == False) & (classification_huvec_filtered['bite'] == False)) |
        ((classification_huvec_filtered['RTS_stage'] == False) & (classification_huvec_filtered['ratio_TSS'] > 1.1)) |
        ((classification_huvec_filtered['bite'] == False) & (
            classification_huvec_filtered['ratio_TSS'] > 1.1))
    )
)
classification_huvec_filtered[
    filter2_1]['structural_category'].value_counts()

# either coverage or all of the other 3
filter2_2 = (
    (classification_huvec_filtered['min_cov'] > 2) | (
        ((classification_huvec_filtered['RTS_stage'] == False) &
         (classification_huvec_filtered['bite'] == False) &
         (classification_huvec_filtered['ratio_TSS'] > 1.1))
    )
)
classification_huvec_filtered[
    filter2_2]['structural_category'].value_counts()


# ISM filter
filter_ism = (
    (classification_huvec_filtered['min_cov'] > 2) | (
        ((classification_huvec_filtered['RTS_stage'] == False) &
         (classification_huvec_filtered['min_cov'] > 0))
    )
)
classification_huvec_filtered[
    filter_ism]['structural_category'].value_counts()


# fusion
filter_fusion = (
    (
        ((classification_huvec_filtered['RTS_stage'] == False) &
         (classification_huvec_filtered['bite'] == False) &
         (classification_huvec_filtered['min_cov'] > 2) &
         (classification_huvec_filtered['ratio_TSS'] > 1.3))
    )
)
classification_huvec_filtered[
    filter_fusion]['structural_category'].value_counts()

# NIC filter
filter_nic = (
    (
        (classification_huvec_filtered['RTS_stage'] == False) &
        (classification_huvec_filtered['bite'] == False) &
        (classification_huvec_filtered['min_cov'] > 0) &
        (classification_huvec_filtered['ratio_TSS'] > 1.1)
    ) | (
        (classification_huvec_filtered['bite'] == False) &
        (classification_huvec_filtered['min_cov'] > 3) &
        (classification_huvec_filtered['ratio_TSS'] > 1.3)
    ) | (
        (classification_huvec_filtered['RTS_stage'] == False) &
        (classification_huvec_filtered['min_cov'] > 4) &
        (classification_huvec_filtered['ratio_TSS'] > 1.3)
    )
)
classification_huvec_filtered[
    filter_nic]['structural_category'].value_counts()

# NNC filter
filter_nnc = (
    (
        (
            (classification_huvec_filtered['RTS_stage'] == False) &
            (classification_huvec_filtered['bite'] == False) &
            (classification_huvec_filtered['min_cov'] > 0) &
            (classification_huvec_filtered['all_canonical'] == 'canonical') &
            (classification_huvec_filtered['ratio_TSS'] > 1.1)
        ) | (
            (classification_huvec_filtered['bite'] == False) &
            (classification_huvec_filtered['min_cov'] > 3) &
            (classification_huvec_filtered['ratio_TSS'] > 1.3)
        ) | (
            (classification_huvec_filtered['RTS_stage'] == False) &
            (classification_huvec_filtered['min_cov'] > 3) &
            (classification_huvec_filtered['all_canonical'] == 'canonical') &
            (classification_huvec_filtered['ratio_TSS'] > 1.3)
        )
    )
)
classification_huvec_filtered[
    filter_nnc]['structural_category'].value_counts()

# genic filter: intron retention
filter_genic = (
    ((classification_huvec_filtered['all_canonical'] == 'canonical') &
     (classification_huvec_filtered['min_cov'] > 2)) |
    (classification_huvec_filtered['min_cov'] > 5)
)
classification_huvec_filtered[
    filter_genic]['structural_category'].value_counts()


# filter rest: intergenic, antisense, genic_intron
filter_rest = (
    (
        (classification_huvec_filtered['all_canonical'] == 'canonical') &
        (classification_huvec_filtered['min_cov'] > 3) &
        (classification_huvec_filtered['RTS_stage'] == False) &
        (classification_huvec_filtered['bite'] == False) &
        (classification_huvec_filtered['ratio_TSS'] > 1.3)
    ) | (
        (classification_huvec_filtered['min_cov'] > 5) &
        (classification_huvec_filtered['RTS_stage'] == False) &
        (classification_huvec_filtered['bite'] == False) &
        (classification_huvec_filtered['ratio_TSS'] > 1.3)
    )
)
classification_huvec_filtered[
    filter_rest]['structural_category'].value_counts()


# ISM filter update
# ratio_TSS quite important: not an degradation artifact, RTS stage not so much

filter_ism_updated = (
    (classification_huvec_filtered['min_cov'] > 2) |
    ((classification_huvec_filtered['ratio_TSS'] > 1.3) &
     (classification_huvec_filtered['min_cov'] > 0))
)
classification_huvec_filtered[
    filter_ism_updated]['structural_category'].value_counts()
