#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Riboseq

UMI_dedup_outdir_transcriptomic=$1

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_13_IP_Puro_1.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_14_IP_Puro_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_13_IP_Puro_1.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_15_IP_Puro_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_15_IP_Puro_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_14_IP_Puro_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_16_IP_CHX_1.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_17_IP_CHX_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_16_IP_CHX_1.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_18_IP_CHX_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_18_IP_CHX_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_17_IP_CHX_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots


python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_10_In_CHX_1.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_11_In_CHX_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots


python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_10_In_CHX_1.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_12_In_CHX_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots


python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_12_In_CHX_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_11_In_CHX_2.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots




python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_07_In_Puro_1.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_08_In_Puro_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_07_In_Puro_1.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_09_In_Puro_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots

python alignments/correlation_plots_transcriptomic.py \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_09_In_Puro_4.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/uf_muellermcnicoll_2025_04_08_In_Puro_3.dedup_idxstats.out \
    ${UMI_dedup_outdir_transcriptomic}/correlation_plots