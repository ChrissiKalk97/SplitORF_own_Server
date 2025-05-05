#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

INDIR="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_1"
OUTDIR_FASTQC1="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/fastqc_unprocessed"
OUTDIR_CUTADAPT="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt"
Umi_adpt_trimmed_path="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/UMI_trimmed_custom"
fastpOut="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim"
fastpFASTQC="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim/fastqc"
Bowtie2_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
Bowtie2_base_name="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/index"
Bowtie2_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"

UMI_Indir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"
UMI_dedup_outdir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/deduplicated"




bash preprocessing_cutadapt_steps.sh \
 $INDIR \
  $OUTDIR_FASTQC1 \
  $OUTDIR_CUTADAPT \
  $Umi_adpt_trimmed_path \
  $fastpOut \
  $fastpFASTQC > "/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/out_reports_of_runs/preprocessing_cutadapt_steps.out" 2>&1

python preprocessing/cutadapt_output_parsing.py \
 "/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/out_reports_of_runs/preprocessing_cutadapt_steps.out" \
  "/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/out_reports_of_runs/cutadapt_summary.csv"



source mapping_steps_ignolia_concatenated_k1.sh \
 $fastpOut \
 $Bowtie2_ref_fasta \
 $Bowtie2_base_name \
 $Bowtie2_out_dir



# deduplicate UMIs
source deduplication/deduplicate_umi_tools.sh \
 $UMI_Indir \
 $UMI_dedup_outdir