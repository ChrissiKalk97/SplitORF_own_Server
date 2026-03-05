#!/bin/bash


# =============================================================================
# Script Name: run_Riboseq_analysis.sh
# Description: This script performs a custom analysis pipeline for the Riboseq 
#               analysis of data created with Vlado's protocol (April 2025):
#              - Step 1: Data preprocessing (skipped)
#              - Step 2: Transcriptomic alignment (Ingolia reference)
#              - Step 3: Transcriptomic deduplication
#              - Step 4: Genomic alignment, deduplication
#              - Step 5: Split-ORF analysis
# Usage:       bash my_pipeline.sh 
# Author:      Christina Kalk
# Date:        2025-05-27
# =============================================================================






eval "$(conda shell.bash hook)"
conda activate Riboseq

################################################################################
# PATH DEFINTIONS                                                              #
################################################################################

genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
ensembl_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
bowtie_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
star_index="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/index"

module_dir="/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_analysis_pipeline"


huvec_dir="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_2_Feb_2026/uf_muellermcnicoll_2026_02_despic_RiboSeq_HUVEC_CM/HUVEC"
cm_dir="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_2_Feb_2026/uf_muellermcnicoll_2026_02_despic_RiboSeq_HUVEC_CM/CM"

out_dir="/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2"
out_dir_cm="${out_dir}"/CM
out_dir_huvec="${out_dir}"/HUVEC

fastp_huvec_dir="${out_dir_huvec}/preprocess/cutadapt/fastp_filter_after_UMI_trim"
fastp_cm_dir="${out_dir_cm}/preprocess/cutadapt/fastp_filter_after_UMI_trim"

report_dir=""${module_dir}"/round_2_huvec_cm_02_2026/outreports_of_runs"
report_dir_cm="${report_dir}"/CM
report_dir_huvec="${report_dir}"/HUVEC

echo "$data_dir"
echo "$out_dir"


######### Preprocessing: cutadapt, UMI trim, fastp #################################
# -d : dedup
# -c : cutadapt
# -p : paired end
# bash "${module_dir}"/Riboseq_pipeline.sh \
#  -d -p -u -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -i "$huvec_dir" -o "$out_dir_huvec" \
#  -r "$report_dir_huvec" -m "${module_dir}"

# bash "${module_dir}"/Riboseq_pipeline.sh \
#  -d -p -u -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -i "$cm_dir" -o "$out_dir_cm" \
#  -r "$report_dir_cm" -m "${module_dir}"

######### Align to transcriptome with bowtie1 #######################################
# bash "${module_dir}"/Riboseq_pipeline.sh \
#  -d -t "${bowtie_ref_fasta}" -i "$huvec_dir" -o "$out_dir_huvec" \
#  -r "$report_dir_huvec" -m "${module_dir}"

# bash "${module_dir}"/Riboseq_pipeline.sh \
#  -d -t "${bowtie_ref_fasta}" -i "$cm_dir" -o "$out_dir_cm" \
#  -r "$report_dir_cm" -m "${module_dir}"


######### Had wrong naming for deduplication bams, correct naming ####################
# bash rename_dedup_bam_files.sh \
# "${out_dir}/CM/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup"

# bash rename_dedup_bam_files.sh \
# "${out_dir}/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup"

# conda activate r-env

# Rscript "${module_dir}"/Ribowaltz/RiboWaltz_Michi_Vlado_2_Ingolia_HUVEC_CM_Feb2026.R \
# "${out_dir}/CM/alignment_concat_transcriptome_Ignolia/Ribowaltz" \
# "${out_dir}/CM/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup"

# Rscript "${module_dir}"/Ribowaltz/RiboWaltz_Michi_Vlado_2_Ingolia_HUVEC_CM_Feb2026.R \
# "${out_dir}/HUVEC/alignment_concat_transcriptome_Ignolia/Ribowaltz" \
# "${out_dir}/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup"


######### Align to genome with STAR and deduplicate #############################
# need to put explicitly the fastp directory, do not need to do this for transcriptomic
# alignment
# bash "${module_dir}"/Riboseq_pipeline.sh \
#  -a "${star_index}" -i "$fastp_huvec_dir" -o "$out_dir_huvec" -g $ensembl_gtf -m "${module_dir}" -q

# bash "${module_dir}"/Riboseq_pipeline.sh \
#  -a "${star_index}" -i "$fastp_cm_dir" -o "$out_dir_cm" -g $ensembl_gtf -m "${module_dir}" -q



# python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_analysis_pipeline/preprocessing/analyze_adapter_dimers_round2_25_02_26.py \
# "${huvec_dir}" "${out_dir_huvec}/preprocess/adapter_dimer_counts_huvec.csv"


######### Creation of metagene plots #############################


# dedup_dir="${out_dir}/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup"
# mkdir -p "${dedup_dir}/coverage_bw_files"
# for bam in ${dedup_dir}/*.bam; do
#     (
#     sample_name=$(basename "$bam" _dedup.bam)

#     bamCoverage \
#     -b "${bam}" \
#     -o "${dedup_dir}/coverage_bw_files/${sample_name}.bw" \
#     --normalizeUsing CPM \
#     --binSize 1
#     )&
# done

# wait

python pyBigWig_metagene_plots.py \
    --path_to_bw_files '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files' \
    --transcript_fai '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai' \
    --mane_transcripts_cds_bed '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed' \
    --condition 'DNOR' \
    --out_path '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files/metagene_plots' \
    --color '#1eb0e6'\
    --region_type 'transcript'

python pyBigWig_metagene_plots.py \
    --path_to_bw_files '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files' \
    --transcript_fai '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai' \
    --mane_transcripts_cds_bed '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed' \
    --condition 'DHYPO' \
    --out_path '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files/metagene_plots' \
    --color '#1eb0e6' \
    --region_type 'transcript'

python pyBigWig_metagene_plots.py \
    --path_to_bw_files '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files' \
    --transcript_fai '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai' \
    --mane_transcripts_cds_bed '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed' \
    --condition 'DNOR' \
    --out_path '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files/metagene_plots' \
    --color '#1eb0e6'\
    --region_type 'cds'

python pyBigWig_metagene_plots.py \
    --path_to_bw_files '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files' \
    --transcript_fai '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai' \
    --mane_transcripts_cds_bed '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed' \
    --condition 'DHYPO' \
    --out_path '/projects/splitorfs/work/Riboseq/Output/HUVEC_CM_round_2/HUVEC/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup/coverage_bw_files/metagene_plots' \
    --color '#1eb0e6' \
    --region_type 'cds'