#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Riboseq

UMI_dedup_outdir_transcriptomic=$1


# correlation plots deduplicated
# control
python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_01_huvec_dnor_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_02_huvec_dnor_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_01_huvec_dnor_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_03_huvec_dnor_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_02_huvec_dnor_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_03_huvec_dnor_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

# hypoxia
python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_04_huvec_dhypo_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_05_huvec_dhypo_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_04_huvec_dhypo_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_06_huvec_dhypo_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_06_huvec_dhypo_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_05_huvec_dhypo_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots




# correlation plots between conditions
python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_01_huvec_dnor_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_04_huvec_dhypo_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots_condition

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_06_huvec_dhypo_4.dedup_idxstats.out  \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_03_huvec_dnor_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots_condition

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_02_huvec_dnor_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_05_huvec_dhypo_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots_condition

# different batches and conditions
python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_01_huvec_dnor_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_05_huvec_dhypo_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots_condition

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_01_huvec_dnor_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_06_huvec_dhypo_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots_condition

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_06_huvec_dhypo_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_03_huvec_dnor_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots_condition