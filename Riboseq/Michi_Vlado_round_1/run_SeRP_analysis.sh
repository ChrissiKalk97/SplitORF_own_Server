#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Riboseq

INDIR="/projects/splitorfs/work/own_data/Riboseq/SeRP/Michi_Vlado_round_1"
OUTDIR_FASTQC1="/projects/splitorfs/work/Riboseq/Output/SeRP/Michi_Vlado_round_1/preprocess/fastqc_unprocessed"
OUTDIR_CUTADAPT="/projects/splitorfs/work/Riboseq/Output/SeRP/Michi_Vlado_round_1/preprocess/cutadapt"
Umi_adpt_trimmed_path="/projects/splitorfs/work/Riboseq/Output/SeRP/Michi_Vlado_round_1/preprocess/cutadapt/UMI_trimmed_custom"
fastpOut="/projects/splitorfs/work/Riboseq/Output/SeRP/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim"
fastpFASTQC="/projects/splitorfs/work/Riboseq/Output/SeRP/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim/fastqc"
Bowtie2_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
Bowtie2_base_name="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/index"
Bowtie2_out_dir="/projects/splitorfs/work/Riboseq/Output/SeRP/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"

UMI_dedup_outdir_transcriptomic=${Bowtie2_out_dir}/deduplicated

UMI_dedup_outdir="/projects/splitorfs/work/Riboseq/Output/SeRP/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"


################################################################################
# QC and PREPROCESSING                                                         #
################################################################################

# bash preprocessing_cutadapt_steps.sh \
#  $INDIR \
#   $OUTDIR_FASTQC1 \
#   $OUTDIR_CUTADAPT \
#   $Umi_adpt_trimmed_path \
#   $fastpOut \
#   $fastpFASTQC > "/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/SeRP/out_reports_of_runs/preprocessing_cutadapt_steps_part2.out" 2>&1

# python preprocessing/cutadapt_output_parsing.py \
#  "/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/SeRP/out_reports_of_runs/preprocessing_cutadapt_steps.out" \
#   "/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/SeRP/out_reports_of_runs/cutadapt_summary.csv"



################################################################################
# TRANSCRIPTOMIC ALIGNMENT                                                     #
################################################################################

# align to transcriptome, still duplicated
# source alignments/bowtie2_align_k1.sh \
#  ${Bowtie2_base_name} \
#  no_index \
#  ${fastpOut} \
#  ${Bowtie2_out_dir} \
#  concat_transcriptome \
#  run_SeRP_analysis_bowtie2.out
# python alignments/analyze_soft_clipping.py ${Bowtie2_out_dir}



################################################################################
# TRANSCRIPTOMIC DEDUPLICATION                                                 #
################################################################################
# deduplicate UMIs
# source deduplication/deduplicate_umi_tools.sh \
#  ${Bowtie2_out_dir} \
#  $UMI_dedup_outdir_transcriptomic \
#  "transcriptomic"



for FQ in "$UMI_dedup_outdir_transcriptomic"/*dedup.bam
do
    fastqc \
    -o "$UMI_dedup_outdir_transcriptomic"/fastqc/ \
    -t 32\
    ${FQ}
done


multiqc --force --filename "$UMI_dedup_outdir_transcriptomic"/fastqc/multiqc_dedup_transcriptomic_alignment "$UMI_dedup_outdir_transcriptomic"/fastqc/



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


# need to change the script
Rscript alignments/PCA_conditions_DeSeq2_SeRP.R $UMI_dedup_outdir_transcriptomic




# align to genome for deduplication
# STAR_index="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/index"
# OutputSTAR="/projects/splitorfs/work/Riboseq/Output/SeRP/Michi_Vlado_round_1/alignment_genome/STAR/only_R1"
# source alignments/genome_alignment_star.sh ${OutputSTAR} ${fastpOut} ${STAR_index}

# python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
#    ${OutputSTAR} \
#     STAR_align_Ribo_genome.csv




# 

# FOR THE DEDUPLICATION: HOW CAN RUN IT IN PARALLEL W/O HAVING TO USE DISTINCT MASTER SCRIPTS?

# # deduplicate UMIs
# source deduplication/deduplicate_umi_tools.sh \
#  $OutputSTAR \
#  $UMI_dedup_outdir \
# genomic



# # filter out secondary and suppl alignments
# files=("${UMI_dedup_outdir}"/*.bam)

# for bam in "${files[@]}"
# do
#     samtools view -F 256 -F 2048 -b ${bam} > \
#      "${$UMI_dedup_outdir}"/$(basename $bam .bam)_filtered.bam
# done