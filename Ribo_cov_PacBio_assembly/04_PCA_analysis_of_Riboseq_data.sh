#!/bin/bash


# =============================================================================
# Script Name: 04_PCA_analysis_of_Riboseq_data.sh
# Description: Perform PCA on each dataset usign featureCounts on genomic mappings:
#              - Set 1: CM - DNOR
#              - Set 2: CM - DHYPO
#              - Step 3: 
#              - Step 4: 
#              - Step 5: 
# Usage:       bash 
# Author:      Christina Wrobel
# Date:        2026-10-09
# =============================================================================



eval "$(conda shell.bash hook)"
conda activate Riboseq



################################################################################
# PATH DEFINTIONS                                                              #
################################################################################
huvec_dir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome"
# dataset with huvec or HUVEC
# dataset BAM in subfolders called endothel
cm_dir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome"
# dataset CM
# dataset ERR in subfolders

huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/HUVEC/HUVEC_LR_SR_support_filtered.gtf"
cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/CM/CM_LR_SR_support_filtered.gtf"


remove_prefix="uf_muellermcnicoll_"
remove_suffix_1=".cutadapt_umi_fastp.only_R1_Aligned.sortedByCoord.out_sorted_dedup_filtered.bam"
remove_suffix_2="_dedup_filtered.bam"
remove_suffix_3="_sorted.bam"
deseq2_design="/home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/04_PCA_analysis_of_Riboseq_data/HUVEC_HYPO_NOR.tsv"

mapping_dir_cm="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/CM_genome"
mapping_dir_huvec="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/PacBio_assemblies_custom_reference_07_10_26/HUVEC_genome"

################################################################################
# GENERATE COUNT FILES FOR ALL SAMPLES TOGETHER FEATURECOUNTS                  #
################################################################################
# for mapping_dir in "$mapping_dir_huvec" "$mapping_dir_cm"; do
#   if [[ "$mapping_dir" == *"CM"* ]];then
#     gtf="$cm_gtf"
#   else
#     gtf="$huvec_gtf"
#   fi

#   shopt -s globstar nullglob
#   bams=( "$mapping_dir"/**/*.bam )
#   echo "Found ${#bams[@]} BAMs:"
#   printf '  %s\n' "${bams[@]}"

#   featureCounts -a "$gtf" -F "GTF" -s 1 --minOverlap 15 -T 32 \
#     -o "$mapping_dir"/gene_mapping_counts.txt "${bams[@]}"

#   awk -v suf="$remove_suffix_1" -v pre="$remove_prefix" 'BEGIN{FS=OFS="\t"}
#     /^Geneid/{ for(i=7;i<=NF;i++){ sub(/.*\//,"",$i); sub(suf"$","",$i); sub("^"pre,"",$i)} }
#     1' "$mapping_dir"/gene_mapping_counts.txt > "$mapping_dir"/gene_mapping_counts.clean.txt

#   awk -v suf="$remove_suffix_2" -v pre="$remove_prefix" 'BEGIN{FS=OFS="\t"}
#     /^Geneid/{ for(i=7;i<=NF;i++){ sub(/.*\//,"",$i); sub(suf"$","",$i); sub("^"pre,"",$i)} }
#     1' "$mapping_dir"/gene_mapping_counts.clean.txt > "$mapping_dir"/gene_mapping_counts.clean2.txt

#   awk -v suf="$remove_suffix_3" -v pre="$remove_prefix" 'BEGIN{FS=OFS="\t"}
#     /^Geneid/{ for(i=7;i<=NF;i++){ sub(/.*\//,"",$i); sub(suf"$","",$i); sub("^"pre,"",$i)} }
#     1' "$mapping_dir"/gene_mapping_counts.clean2.txt > "$mapping_dir"/gene_mapping_counts.clean3.txt
# done

conda activate r-env
Rscript /home/ckalk/scripts/SeRP/PCA_conditions_CDS_generalized.R \
 "$mapping_dir_cm"/gene_mapping_counts.clean3.txt \
 /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/04_PCA_analysis_of_Riboseq_data/CM_HYPO_NOR.tsv

 Rscript /home/ckalk/scripts/SeRP/PCA_conditions_CDS_generalized.R \
 "$mapping_dir_huvec"/gene_mapping_counts.clean3.txt \
 /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/04_PCA_analysis_of_Riboseq_data/HUVEC_HYPO_NOR.tsv

 Rscript /home/ckalk/scripts/SeRP/PCA_conditions_CDS_generalized.R \
 "$mapping_dir_huvec"/gene_mapping_counts.clean3.txt \
 /home/ckalk/scripts/SplitORFs/Ribo_cov_PacBio_assembly/04_PCA_analysis_of_Riboseq_data/endothel_IL1beta_control.tsv