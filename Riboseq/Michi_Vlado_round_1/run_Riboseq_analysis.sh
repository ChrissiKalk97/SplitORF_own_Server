#!/bin/bash


# =============================================================================
# Script Name: run_Riboseq_analysis.sh
# Description: This script performs a custom analysis pipeline for the Riboseq 
#               analysis of data created with Vlado's protocol (April 2025):
#              - Step 1: Data preprocessing
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

INDIR="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_1"
OUTDIR_FASTQC1="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/fastqc_unprocessed"
OUTDIR_CUTADAPT="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt"
Umi_adpt_trimmed_path="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/UMI_trimmed_custom"
fastpOut="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim"
fastpFASTQC="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim/fastqc"
Bowtie2_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
Bowtie2_base_name="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/index"
Bowtie2_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"


UMI_Indir_transcriptomic="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"
UMI_dedup_outdir_transcriptomic="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/deduplicated"


################################################################################
# QC and PREPROCESSING                                                         #
################################################################################

# bash preprocessing_cutadapt_steps.sh \
#  $INDIR \
#   $OUTDIR_FASTQC1 \
#   $OUTDIR_CUTADAPT \
#   $Umi_adpt_trimmed_path \
#   $fastpOut \
#   $fastpFASTQC > "/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/out_reports_of_runs/preprocessing_cutadapt_steps.out" 2>&1

# python preprocessing/cutadapt_output_parsing.py \
#  "/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/out_reports_of_runs/preprocessing_cutadapt_steps.out" \
#   "/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/out_reports_of_runs/cutadapt_summary.csv"

################################################################################
# TRANSCRIPTOMIC ALIGNMENT                                                     #
################################################################################

# source alignments/bowtie2_align_k1_only_R1.sh \
#  ${Bowtie2_base_name} \
#  no_index \
#  ${fastpOut} \
#  ${Bowtie2_out_dir} \
#  concat_transcriptome \
#  transcriptomic_mapping_k1_R1.out \
#  > transcriptomic_mapping_k1_R1.out 2>&1


# # count soft clipping in transcriptomic alignments
# python alignments/analyze_soft_clipping.py ${Bowtie2_out_dir}



################################################################################
# TRANSCRIPTOMIC DEDUPLICATION                                                 #
################################################################################
# deduplicate UMIs
# source deduplication/deduplicate_umi_tools.sh \
#  ${UMI_Indir_transcriptomic}/filtered/q10 \
#  ${UMI_Indir_transcriptomic}/filtered/q10/dedup \
#  transcriptomic


# source deduplication/deduplicate_umi_tools.sh \
#  ${UMI_Indir_transcriptomic}/filtered \
#  $UMI_dedup_outdir_transcriptomic/ \
#  transcriptomic



# count soft clipping in transcriptomic alignments
# python alignments/analyze_soft_clipping.py $UMI_dedup_outdir_transcriptomic

if [ ! -d ${UMI_Indir_transcriptomic}/filtered/q10/dedup/DEGs ]; then
    mkdir ${UMI_Indir_transcriptomic}/filtered/q10/dedup/DEGs
fi
# Rscript alignments/PCA_conditions_DeSeq2.R ${UMI_Indir_transcriptomic}/filtered/q10/dedup

# conda activate pygtftk
# # obtain gene IDs of differentially expressed transcripts
# python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
#  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
#  ${UMI_Indir_transcriptomic}/filtered/q10/dedup/DEGs/hypo_vs_norm_Riboseq_DEGs_1.txt \
#  ${UMI_Indir_transcriptomic}/filtered/q10/dedup/DEGs/hypo_vs_norm_Riboseq_DEGs_gene_ids_1.txt

#  python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
#  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
#  ${UMI_Indir_transcriptomic}/filtered/q10/dedup/DEGs/hypo_vs_norm_downreg_Riboseq_DEGs_1.txt \
#  ${UMI_Indir_transcriptomic}/filtered/q10/dedup/DEGs/hypo_vs_norm_downreg_Riboseq_DEGs_gene_ids_1.txt

#   python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
#  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
#  ${UMI_Indir_transcriptomic}/filtered/q10/dedup/DEGs/hypo_vs_norm_up_down_Riboseq_DEGs_1.txt \
#  ${UMI_Indir_transcriptomic}/filtered/q10/dedup/DEGs/hypo_vs_norm_up_down_Riboseq_DEGs_gene_ids_1.txt


# # FILTER TO KEEP ONLY MRNA MAPPING READS FOR DOWNSTREAM PURPOSES
# # MOVE THE SCATTERPLOTS HERE
# bash alignments/create_correlation_plots_Riboseq.sh ${UMI_Indir_transcriptomic}/filtered/q10/dedup

if [ ! -d ${UMI_Indir_transcriptomic}/Ribowaltz ]; then
    mkdir ${UMI_Indir_transcriptomic}/Ribowaltz
fi

# Rscript Ribowaltz/RiboWaltz_Michi_Vlado_1_Ignolia_single_samples.R 

################################################################################
# GENOMIC ALIGNMENT                                                            #
################################################################################

STAR_index="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/index"
OutputSTAR="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1"
# source alignments/genome_alignment_star.sh ${OutputSTAR} ${fastpOut} ${STAR_index}

Ens110_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_rRNA/Ens110"
NCBI_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_rRNA/NCBI"

# source alignments/bowtie2_align_k1_only_R1.sh \
#  ${Ens110_dir}/index \
#  /projects/splitorfs/work/Riboseq/data/contamination/rRNA_transcripts_110.fasta \
#  ${fastpOut} \
#  $Ens110_dir \
#  EnsrRNA \
#  None

#  source alignments/bowtie2_align_k1_only_R1.sh \
#  ${NCBI_dir}/index \
#  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/NCBI_rRNAs.fasta \
#  ${fastpOut} \
#  $NCBI_dir \
#  NCBIrRNA \
#  None
 


# python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
#     /projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1 \
#     STAR_align_Ribo_genome.csv


# python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
#     /projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR \
#     STAR_align_Ribo_genome.csv



UMI_Indir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1"
UMI_dedup_outdir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"


# deduplicate UMIs
# source deduplication/deduplicate_umi_tools.sh \
#  $UMI_Indir \
#  $UMI_dedup_outdir



# filter out secondary and suppl alignments
files=("${UMI_dedup_outdir}"/*_dedup.bam)

for bam in "${files[@]}"
do
    samtools view -F 256 -F 2048 -q 10 -b ${bam} > \
     "${UMI_dedup_outdir}"/$(basename $bam .bam)_filtered.bam

     samtools index "${UMI_dedup_outdir}"/$(basename $bam .bam)_filtered.bam

     samtools idxstats "${UMI_dedup_outdir}"/$(basename $bam .bam)_filtered.bam > \
    "${UMI_dedup_outdir}"/$(basename $bam .bam)_filtered_idxstats.out

    samtools stats "${UMI_dedup_outdir}"/$(basename $bam .bam)_filtered.bam > \
    "${UMI_dedup_outdir}"/$(basename $bam .bam)_filtered_stats.out

    samtools flagstat "${UMI_dedup_outdir}"/$(basename $bam .bam)_filtered.bam > \
    "${UMI_dedup_outdir}"/$(basename $bam .bam)_filtered_flagstat.out
done


# # # analyze soft clipping of genomic deduplciated reads
# python alignments/analyze_soft_clipping.py $UMI_dedup_outdir


# run FeatureCounts to get mapping percentages
Rscript preprocessing/analyze_mappings/genome_aligned_reads_biotype_counting.R




################################################################################
# SO ANALYSIS OF GENOMIC ALIGNMENTS                                            #
################################################################################
bash Riboseq_SO_empirical_intersection/wrapper_empirical_intersection_random_resample.sh