#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Riboseq

preprocess_script_dir="/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_analysis_pipeline"

indir="/projects/serp/work/data/SeRP_March_2026/Importins"
outdir="/projects/serp/work/Output/March_2026/Importins"
outdir_fastqc="${outdir}/preprocess/fastqc_unprocessed"
outdir_cutadapt="${outdir}/preprocess/cutadapt"
fastp_out="${outdir}/preprocess/fastp"
fastp_fastqc="${outdir}/preprocess/fastp/fastqc"
transcriptome_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
Genome_Fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
EnsemblFilteredRef="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"

Bowtie1_base_name="${outdir}/transcriptome_mapping/bowtie1/index"
Bowtie1_out_dir="${outdir}/transcriptome_mapping/bowtie1"

report_dir="/home/ckalk/scripts/SplitORFs/Riboseq/SeRP/Importins_HA_March_2026/outreports_of_runs"

analysis_script_dir="/home/ckalk/scripts/SplitORFs/Riboseq/SeRP"
coverage_script_dir="${analysis_script_dir}/coverage_plots"
mane_gtf="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf"

mkdir -p "/projects/serp/work/Output/March_2026"
mkdir -p "${outdir}"
mkdir -p "${outdir}/preprocess"
mkdir -p "${outdir}/transcriptome_mapping"

################################################################################
# QC and PREPROCESSING                                                         #
################################################################################

# bash "${analysis_script_dir}"/preprocessing_cutadapt_steps_importins.sh \
#  "$indir" \
#   "$outdir_fastqc" \
#   "$outdir_cutadapt" \
#   "$fastp_out" \
#   "$fastp_fastqc" \
#  "$preprocess_script_dir" \
# > "${report_dir}/preprocessing_cutadapt_importins.out" 2>&1

# python "${preprocess_script_dir}/preprocessing/cutadapt_output_parsing.py" \
#  "${report_dir}/preprocessing_cutadapt_importins.out" \
#   "${report_dir}/cutadapt_summary.csv"



################################################################################
# TRANSCRIPTOMIC ALIGNMENT                                                     #
################################################################################
source "${preprocess_script_dir}"/alignments/bowtie1_align_21_10_25.sh \
 "$Bowtie1_base_name" \
 no_index \
 "$fastp_out" \
 "$Bowtie1_out_dir" \
 concat_transcriptome \
 "${report_dir}"/run_SeRP_analysis_bowtie1.out \
 "$preprocess_script_dir" \
 "$indir" #\
#  > "${report_dir}"/run_SeRP_analysis_bowtie1.out 2>&1
# "$transcriptome_fasta" \


# if [ ! -d ${Bowtie1_out_dir}/filtered/q10/DEGs ]; then
#     mkdir ${Bowtie1_out_dir}/filtered/q10/DEGs
# fi

# conda activate r-env

# Rscript  "${analysis_script_dir}"/Importins_HA_March_2026/PCA_conditions_DeSeq2_SeRP_importins_HA.R \
# "${Bowtie1_out_dir}"/filtered/q10


# bash map_DEG_tID_to_gID.sh "${Bowtie1_out_dir}"


# conda activate Riboseq

# if [ ! -d ${Bowtie1_out_dir}/filtered/q10/Ribowaltz ]; then
#     mkdir ${Bowtie1_out_dir}/filtered/q10/Ribowaltz
# fi

# Rscript RiboWaltz_SeRP_importins_single_samples_bowtie1.R \
# "${Bowtie1_out_dir}"/filtered/q10



# ################################################################################
# # SeRP coverage plots                                                          #
# ################################################################################
# if [[ ! -d "${Bowtie1_out_dir}"/filtered/q10/enrichment_plots_CDS ]]; then
#     mkdir "${Bowtie1_out_dir}"/filtered/q10/enrichment_plots_CDS
# fi 

# bash ${coverage_script_dir}/create_coverage_plots_codons_CDS_different_buffer_transcript.sh \
#     "${Bowtie1_out_dir}"/filtered/q10 \
#     ""${Bowtie1_out_dir}"/filtered/q10/enrichment_plots_CDS/CDS_coordinates" \
#     ${coverage_script_dir} \
#     $mane_gtf
















