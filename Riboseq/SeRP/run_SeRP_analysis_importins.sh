#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Riboseq

script_dir="/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1"

INDIR="/projects/serp/work/data/SeRP_April_2025/importins"
OUTDIR_FASTQC1="/projects/serp/work/Output/April_2025/importins/preprocess/fastqc_unprocessed"
OUTDIR_CUTADAPT="/projects/serp/work/Output/April_2025/importins/preprocess/cutadapt"
fastpOut="/projects/serp/work/Output/April_2025/importins/preprocess/fastp"
fastpFASTQC="/projects/serp/work/Output/April_2025/importins/preprocess/fastp/fastqc"
Bowtie2_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
Bowtie2_base_name="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/index"
Bowtie2_out_dir="/projects/serp/work/Output/April_2025/importins/transcriptome_mapping"
Genome_Fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
EnsemblFilteredRef="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"

Bowtie1_base_name="/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/bowtie1/index"
Bowtie1_out_dir="/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/bowtie1"

coverage_script_dir="/home/ckalk/scripts/SplitORFs/Riboseq/SeRP/coverage_plots"
mane_gtf="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf"

################################################################################
# QC and PREPROCESSING                                                         #
################################################################################

# bash preprocessing_cutadapt_steps_importins.sh \
#  $INDIR \
#   $OUTDIR_FASTQC1 \
#   $OUTDIR_CUTADAPT \
#   $fastpOut \
#   $fastpFASTQC\
# > "out_reports_of_runs/preprocessing_cutadapt_importins.out" 2>&1

# python ${script_dir}/preprocessing/cutadapt_output_parsing.py \
#  "out_reports_of_runs/preprocessing_cutadapt_importins.out" \
#   "/home/ckalk/scripts/SplitORFs/Riboseq/SeRP/out_reports_of_runs/cutadapt_summary.csv"



################################################################################
# TRANSCRIPTOMIC ALIGNMENT                                                     #
################################################################################
# align to transcriptome
# source ${script_dir}/alignments/bowtie2_align_k1_only_R1.sh \
#  ${Bowtie2_base_name} \
#  index \
#  ${fastpOut} \
#  ${Bowtie2_out_dir} \
#  concat_transcriptome \
#  /home/ckalk/scripts/SplitORFs/Riboseq/SeRP/out_reports_of_runs/run_SeRP_analysis_bowtie2.out \
#  > /home/ckalk/scripts/SplitORFs/Riboseq/SeRP/out_reports_of_runs/run_SeRP_analysis_bowtie2.out 2>&1

# source ${script_dir}/alignments/bowtie1_align_21_10_25.sh \
#  ${Bowtie1_base_name} \
#  ${Bowtie2_ref_fasta} \
#  ${fastpOut} \
#  ${Bowtie1_out_dir} \
#  concat_transcriptome \
#  /home/ckalk/scripts/SplitORFs/Riboseq/SeRP/out_reports_of_runs/run_SeRP_analysis_bowtie1.out \
#  > /home/ckalk/scripts/SplitORFs/Riboseq/SeRP/out_reports_of_runs/run_SeRP_analysis_bowtie1.out 2>&1




# python ${script_dir}/alignments/analyze_soft_clipping.py ${Bowtie1_out_dir}/filtered



# # bash create_correlation_plots_reps_SeRP_importins.sh \
# #  "${Bowtie1_out_dir}"/filtered/q10

if [ ! -d ${Bowtie1_out_dir}/filtered/q10/DEGs ]; then
    mkdir ${Bowtie1_out_dir}/filtered/q10/DEGs
fi

Rscript  PCA_conditions_DeSeq2_SeRP_importins.R \
"${Bowtie1_out_dir}"/filtered/q10


bash map_DEG_tID_to_gID.sh "${Bowtie1_out_dir}"


conda activate Riboseq

if [ ! -d ${Bowtie1_out_dir}/filtered/q10/Ribowaltz ]; then
    mkdir ${Bowtie1_out_dir}/filtered/q10/Ribowaltz
fi

Rscript RiboWaltz_SeRP_importins_single_samples_bowtie1.R \
"${Bowtie1_out_dir}"/filtered/q10



# ################################################################################
# # SeRP coverage plots                                                          #
# ################################################################################
# # bash create_coverage_plots_importins_over_mock.sh \
# #  "${Bowtie2_out_dir}"/filtered/q10

# bash ${coverage_script_dir}/create_coverage_plots_codons_CDS_different_buffer_transcript.sh \
#     "${Bowtie1_out_dir}"/filtered/q10 \
#     ""${Bowtie1_out_dir}"/filtered/q10/enrichment_plots_CDS/CDS_coordinates" \
#     ${coverage_script_dir} \
#     $mane_gtf

################################################################################
# Check al counts with bam multicov                                            #
################################################################################

# faidx --transform bed ${Bowtie2_ref_fasta} > $(dirname ${Bowtie2_ref_fasta})/Ignolia_transcriptome_and_contamination.bed

# reference_bed=$(dirname ${Bowtie2_ref_fasta})/Ignolia_transcriptome_and_contamination.bed

# for bam in "${Bowtie2_out_dir}"/filtered/q10/*.bam; do
    # bedtools multicov -bams $bam -bed $reference_bed > "${Bowtie2_out_dir}"/filtered/q10/$(basename $bam .bam)_mulitcov.bed
#     diff <(cut -f 1,4 "${Bowtie2_out_dir}"/filtered/q10/$(basename $bam .bam)_mulitcov.bed)\
#      <(cut -f 1,3  "${Bowtie2_out_dir}"/filtered/q10/$(basename $bam .bam)_idxstats.out)
# done














