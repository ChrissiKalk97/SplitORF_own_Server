#!/bin/bash

#----- This script runs gffcompare on 2 assemblies to merge and compare them ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio

reference_gtf="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"

mando_rescued_cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf"
stringtie_cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3/CM/CM_strigntie3_assembly_renamed_filtered.gtf"
outdir_cm="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/CM"
prefix_cm="CM_mando_stringtie_combined"

mando_rescued_huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf"
stringtie_huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3/HUVEC/HUVEC_strigntie3_assembly_renamed_filtered.gtf"
outdir_huvec="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/HUVEC"
prefix_huvec="HUVEC_mando_stringtie_combined"

genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"


#################################################################################
# ------------------ MERGE WITH GFFCOMPARE                   ------------------ #
#################################################################################

# bash gffcompare_stringtie_mando.sh \
#     $reference_gtf \
#     $mando_rescued_cm_gtf \
#     $stringtie_cm_gtf \
#     $outdir_cm \
#     $prefix_cm


# bash gffcompare_stringtie_mando.sh \
#     $reference_gtf \
#     $mando_rescued_huvec_gtf \
#     $stringtie_huvec_gtf \
#     $outdir_huvec \
#     $prefix_huvec



# script_dir="/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3"
# genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_cm/CM_mando_stringtie_combined.combined.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_cm/SQANTI_QC

#  bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_huvec/HUVEC_mando_stringtie_combined.combined.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_huvec/SQANTI_QC


#################################################################################
# ------------------ RUN SPLIT-ORFs PIPELINE                 ------------------ #
#################################################################################
bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
 $outdir_tama/CM/CM_merged_tama_gene_id.gtf \
 /home/ckalk/tools/SplitORF_pipeline \
 $genome_fasta


bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
 $outdir_tama/HUVEC/HUVEC_merged_tama_gene_id.gtf \
 /home/ckalk/tools/SplitORF_pipeline \
 $genome_fasta

#################################################################################
# ------------------ RUN FIFTYNT PIPELINE                    ------------------ #
#################################################################################
# bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_fiftynt_on_assembly.sh \
#     $out_path/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf \
#     /home/ckalk/tools/NMD_fetaure_composition \
#     $genome_fasta \
#     $ensembl_full_gtf \
#     CM_mando_rescued_50nt.csv


# bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_fiftynt_on_assembly.sh \
#     $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf \
#     /home/ckalk/tools/NMD_fetaure_composition \
#     $genome_fasta \
#     $ensembl_full_gtf \
#     HUVEC_mando_rescued_50nt.csv
