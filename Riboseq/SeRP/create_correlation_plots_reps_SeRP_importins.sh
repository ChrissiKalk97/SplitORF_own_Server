#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Riboseq

bowtie2_alns_q10_filtered_dir=$1

if [ ! -d ${bowtie2_alns_q10_filtered_dir}/correlation_plots ]; then
        mkdir ${bowtie2_alns_q10_filtered_dir}/correlation_plots
fi

# CHX A2
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_18_RR_A2_CHX_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_19_RR_A2_CHX_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_18_RR_A2_CHX_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_20_RR_A2_CHX_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_18_RR_A2_CHX_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_21_RR_A2_CHX_E4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_19_RR_A2_CHX_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_20_RR_A2_CHX_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_19_RR_A2_CHX_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_21_RR_A2_CHX_E4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_20_RR_A2_CHX_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_21_RR_A2_CHX_E4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots


# CHX B1
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_25_RR_B1_CHX_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_26_RR_B1_CHX_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_25_RR_B1_CHX_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_27_RR_B1_CHX_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_25_RR_B1_CHX_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_28_RR_B1_CHX_E4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_26_RR_B1_CHX_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_27_RR_B1_CHX_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_26_RR_B1_CHX_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_28_RR_B1_CHX_E4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_27_RR_B1_CHX_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_28_RR_B1_CHX_E4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots


# CHX In
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_32_RR_In_CHX_1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_33_RR_In_CHX_2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_32_RR_In_CHX_1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_34_RR_In_CHX_3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_32_RR_In_CHX_1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_35_RR_In_CHX_4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_33_RR_In_CHX_2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_34_RR_In_CHX_3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_33_RR_In_CHX_2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_35_RR_In_CHX_4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_34_RR_In_CHX_3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_35_RR_In_CHX_4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

# Puro A2
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_22_RR_A2_Puro_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_23_RR_A2_Puro_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_22_RR_A2_Puro_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_24_RR_A2_Puro_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots


python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_23_RR_A2_Puro_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_24_RR_A2_Puro_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots


# Puro B1
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_29_RR_B1_Puro_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_30_RR_B1_Puro_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_29_RR_B1_Puro_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_31_RR_B1_Puro_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots


python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_30_RR_B1_Puro_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_31_RR_B1_Puro_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

# Puro In
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_36_RR_In_Puro_1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_37_RR_In_Puro_2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_36_RR_In_Puro_1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_38_RR_In_Puro_3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots


python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_37_RR_In_Puro_2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_38_RR_In_Puro_3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots


# CHX Mock
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_39_RR_M_CHX_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_40_RR_M_CHX_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_39_RR_M_CHX_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_41_RR_M_CHX_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_39_RR_M_CHX_E1.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_42_RR_M_CHX_E4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_40_RR_M_CHX_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_41_RR_M_CHX_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_40_RR_M_CHX_E2.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_42_RR_M_CHX_E4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/correlation_plots_transcriptomic.py \
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_41_RR_M_CHX_E3.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/uf_muellermcnicoll_2025_05_42_RR_M_CHX_E4.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10_idxstats.out\
    ${bowtie2_alns_q10_filtered_dir}/correlation_plots
