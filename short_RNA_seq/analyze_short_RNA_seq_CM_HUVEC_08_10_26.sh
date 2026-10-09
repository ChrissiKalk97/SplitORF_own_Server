#!/bin/bash

#----- This script performs the preprocessing steps for the short RNA-seq samples ----- #
# ----- first cutadapt is used to trim adapters and filter for quality            ----- #
# ----- FASTQC is run to obtain quality control metrics and plots                 ----- #

eval "$(conda shell.bash hook)"
conda activate Riboseq

raw_data_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData"
merged_data_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/merged"
raw_data_fastqc_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData/fastqc"

multiQC_outname="multiqc"
# 40 FASTQ files are indicated in the MD5 checksums

outidr_fastp="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_fastp"
Genome_Fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
Decoys="/projects/splitorfs/work/reference_files/decoys.txt"

out_dir_huvec="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_TAMA_08_10_26"
STAR_out_dir_huvec="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_TAMA_08_10_26/STAR"
salmon_out_dir_huvec_tama="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Salmon/HUVEC_TAMA_08_10_26"
huvec_tama_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/HUVEC/HUVEC_LR_SR_support_filtered.gtf"

mkdir -p "/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Salmon"
mkdir -p "$out_dir_huvec"

# bash fastqc_multiqc_RNA_seq.sh ${raw_data_dir_huvec} ${raw_data_fastqc_dir_huvec} ${multiQC_outname} raw

# bash merge_fastq_files.sh ${raw_data_dir_huvec} ${merged_data_dir_huvec}

# bash fastqc_multiqc_RNA_seq.sh ${merged_data_dir_huvec} ${merged_data_dir_huvec}/fastqc ${multiQC_outname} raw

# bash trim_adapters_RNA_seq.sh ${merged_data_dir_huvec} ${outidr_fastp} _fastp

# bash fastqc_multiqc_RNA_seq.sh ${outidr_fastp} ${outidr_fastp}/fastqc ${multiQC_outname} fastp

bash STAR_map_short_RNA_seq.sh \
 $outidr_fastp \
 $STAR_out_dir_huvec \
 $STAR_out_dir_huvec/index \
 $Genome_Fasta \
 ${huvec_tama_gtf}


python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
    /projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/STAR \
    STAR_align_genome.csv


bash Salmon_quantification.sh \
  ${Genome_Fasta} \
  ${huvec_tama_gtf} \
  ${salmon_out_dir_huvec_tama} \
  ${outidr_fastp} \
  ${salmon_out_dir_huvec_tama} \
  ${Decoys}



#################################################################################
# ------------------ ANALYZE CM SHORT-READ DATA              ------------------ #
#################################################################################
raw_data_dir_cm="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/CM_short_reads/X208SC25032333-Z01-F002/01.RawData"
merged_data_dir_cm="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/CM_short_reads/merged"
raw_data_fastqc_dir_cm="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/CM_short_reads/X208SC25032333-Z01-F002/01.RawData/fastqc"

multiQC_outname="multiqc"
# 40 FASTQ files are indicated in the MD5 checksums

out_dir_cm="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/CM_TAMA_08_10_26"
STAR_out_dir_cm="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/CM_TAMA_08_10_26/STAR"
salmon_out_dir_cm_tama="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Salmon/CM_TAMA_08_10_26"
Genome_Fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"


cm_tama_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/CM/CM_LR_SR_support_filtered.gtf"


mkdir -p "$out_dir_cm"


# bash fastqc_multiqc_RNA_seq.sh ${raw_data_dir_cm} ${raw_data_fastqc_dir_cm} ${multiQC_outname} raw

# bash merge_fastq_files.sh ${raw_data_dir_cm} ${merged_data_dir_cm}

# bash fastqc_multiqc_RNA_seq.sh ${merged_data_dir_cm} ${merged_data_dir_cm}/fastqc ${multiQC_outname} raw

bash trim_adapters_RNA_seq.sh ${merged_data_dir_cm} ${outidr_fastp_cm} _fastp

bash fastqc_multiqc_RNA_seq.sh ${outidr_fastp} ${outidr_fastp_cm}/fastqc ${multiQC_outname} fastp

bash STAR_map_short_RNA_seq.sh \
 $outidr_fastp_cm \
 $STAR_out_dir_cm \
 $STAR_out_dir_cm/index \
 $Genome_Fasta \
 ${cm_tama_gtf}


python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
    /projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/STAR \
    STAR_align_genome.csv


bash Salmon_quantification.sh \
  ${Genome_Fasta} \
  ${cm_tama_gtf} \
  ${salmon_out_dir_cm_tama} \
  ${outidr_fastp_cm} \
  ${salmon_out_dir_cm_tama} \
  ${Decoys}