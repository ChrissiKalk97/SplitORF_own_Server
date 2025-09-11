#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

# ensembl_filtered_gtf=$1
# genome_fasta=$2
# consensus_reads_fofn=$3


ensembl_filtered_gtf="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"
ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
consensus_reads_fofn_HUVEC="pacbio_consensus_HUVEC.fofn"
consensus_reads_fofn_CM="./pacbio_consensus_CM.fofn"
out_path="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters"
# with the default filter settings
out_path_filter="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion"
bam_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"


if [[ ! -d "$out_path" ]]; then
    mkdir $out_path
fi

if [[ ! -d "$out_path/HUVEC" ]]; then
    mkdir $out_path/HUVEC
fi

if [[ ! -d "$out_path/CM" ]]; then
    mkdir $out_path/CM
fi

if [[ ! -d "$bam_dir/fastq" ]]; then
    mkdir $bam_dir/fastq
fi




shopt -s nullglob
bam_files=("${bam_dir}"/*bam)


#################################################################################
# ------------------ GET FASTQ FILES                         ------------------ #
#################################################################################
# for bam in "${bam_files[@]}"; 
# do
#     sample=$(basename $bam .bam)
#     bam2fastq -u -o $bam_dir/fastq/$sample $bam
# done

#################################################################################
# ------------------ RUN MANDOLORION TO CREATE ASSEMBLY      ------------------ #
#################################################################################
# python3 /home/ckalk/scripts/SplitORFs/PacBio_analysis/tools/Mandalorion/Mando.py \
# -p $out_path/HUVEC \
# -t 32 \
# -g $ensembl_filtered_gtf \
# -G $genome_fasta \
# --minimum_ratio 0 \
# --minimum_reads 2 \
# --minimum_feature_count 2 \
# --Acutoff 1 \
# -f /projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_50NMD_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_5NMD_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_DHYPO_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_DMSO_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_DNOR_merged_lima_refined.fastq

# python3 /home/ckalk/scripts/SplitORFs/PacBio_analysis/tools/Mandalorion/Mando.py \
# -p $out_path/CM \
# -t 32 \
# -g $ensembl_filtered_gtf \
# -G $genome_fasta \
# --minimum_ratio 0 \
# --minimum_reads 2 \
# --minimum_feature_count 2 \
# --Acutoff 1 \
# -f /projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/CM_5NMD_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/CM_DHYPO_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/CM_DNOR_merged_lima_refined.fastq

# get bam files of minimap2 alignments made with Mando
# samtools view -bo $out_path/HUVEC/tmp/Isoforms.aligned.out.filtered.bam \
#  $out_path/HUVEC/tmp/Isoforms.aligned.out.filtered.sam 
# samtools view -bo $out_path/CM/tmp/Isoforms.aligned.out.filtered.bam \
#  $out_path/CM/tmp/Isoforms.aligned.out.filtered.sam 

# samtools sort $out_path/HUVEC/tmp/Isoforms.aligned.out.filtered.bam \
#  -o $out_path/HUVEC/tmp/Isoforms.aligned.out.filtered.sorted.bam
# samtools sort $out_path/CM/tmp/Isoforms.aligned.out.filtered.bam \
#  -o $out_path/CM/tmp/Isoforms.aligned.out.filtered.sorted.bam

# samtools index $out_path/HUVEC/tmp/Isoforms.aligned.out.filtered.sorted.bam
# samtools index $out_path/CM/tmp/Isoforms.aligned.out.filtered.sorted.bam


# #################################################################################
# # ------------RENAME ASSEMBLY TO GET GENE ID AND TID ENSEMBL ------------------ #
# #################################################################################
# python mandalorion/rename_gene_id_name_mando_gtf.py $out_path/HUVEC/Isoforms.filtered.clean.gtf $out_path/HUVEC/HUVEC_mando_gene_id.gtf
# python mandalorion/rename_gene_id_name_mando_gtf.py $out_path/CM/Isoforms.filtered.clean.gtf $out_path/CM/CM_mando_gene_id.gtf



# gffread $out_path/HUVEC/HUVEC_mando_gene_id.gtf \
#  -g $genome_fasta -w $out_path/HUVEC/HUVEC_mando_gene_id.fasta

# gffread $out_path/CM/CM_mando_gene_id.gtf \
#  -g $genome_fasta -w $out_path/CM/CM_mando_gene_id.fasta


#################################################################################
# ------------------ LR support/expression of isoforms       ------------------ #
#################################################################################
# python mandalorion/plot_isoform_quantification_mando.py $out_path/HUVEC 5 50
# python mandalorion/plot_isoform_quantification_mando.py $out_path/CM 3 50
# python mandalorion/plot_isoform_quantification_mando.py $out_path_filter/HUVEC 5 50
# python mandalorion/plot_isoform_quantification_mando.py $out_path_filter/CM 3 50




#################################################################################
# ------------------ fl counts for SQANTI3                   ------------------ #
#################################################################################
# conda activate Riboseq
# python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/sqanti3/get_fl_count_from_mando_quant_output.py \
#   $out_path/HUVEC

# python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/sqanti3/get_fl_count_from_mando_quant_output.py \
#   $out_path/CM



#################################################################################
# ------------------ Run SQANTI                              ------------------ #
#################################################################################
# conda activate pacbio
# bash mandalorion/assess_mando_sqanti3/run_sqanti3_on_mando_updated_arguments_august_2025.sh \
#  $genome_fasta \
#  $ensembl_filtered_gtf \
#  $out_path \
#  "/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Mandalorion_raw_updated_parameters" \
#  "/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/fastp" \
#  "/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3"


#################################################################################
# ------------------ RUN SPLIT-ORFs PIPELINE                 ------------------ #
#################################################################################
# bash SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
# $out_path/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf \
# /home/ckalk/tools/SplitORF_pipeline \
# $genome_fasta


# bash SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
# $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf \
# /home/ckalk/tools/SplitORF_pipeline \
# $genome_fasta

#################################################################################
# ------------------ RUN FIFTYNT PIPELINE                    ------------------ #
#################################################################################
bash SplitORF_scripts/run_fiftynt_on_assembly.sh \
    $out_path/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf \
    /home/ckalk/tools/NMD_fetaure_composition \
    $genome_fasta \
    $ensembl_full_gtf \
    CM_mando_rescued_50nt.csv


bash SplitORF_scripts/run_fiftynt_on_assembly.sh \
    $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf \
    /home/ckalk/tools/NMD_fetaure_composition \
    $genome_fasta \
    $ensembl_full_gtf \
    HUVEC_mando_rescued_50nt.csv



#################################################################################
# ------------------ COMPARE TO ENSEMBL FILTERED ASSEMBLY    ------------------ #
#################################################################################

# if [[ ! -d "$out_path/HUVEC/gffcompare_only_r" ]]; then
#     mkdir $out_path/HUVEC/gffcompare_only_r
# fi

# gffcompare -o $out_path/HUVEC/gffcompare_only_r/HUVEC_gffcompare_only_r\
#  -r $ensembl_filtered_gtf\
#   $out_path/HUVEC/HUVEC_mando_gene_id.gtf


# if [[ ! -d "$out_path/CM/gffcompare_only_r" ]]; then
#     mkdir $out_path/CM/gffcompare_only_r
# fi

# gffcompare -o $out_path/CM/gffcompare_only_r/CM_gffcompare_only_r\
#  -r $ensembl_filtered_gtf\
#   $out_path/CM/CM_mando_gene_id.gtf


# python mandalorion/get_equal_ejc_isoforms.py \
#  $out_path/CM/gffcompare_only_r/CM_gffcompare_only_r.CM_mando_gene_id.gtf.tmap

# python mandalorion/get_equal_ejc_isoforms.py \
#  $out_path/HUVEC/gffcompare_only_r/HUVEC_gffcompare_only_r.HUVEC_mando_gene_id.gtf.tmap

#  python mandalorion/count_nr_novel_nmd_transcripts.py \
#  /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/CM/fifty_nt_result/CM_mando_gene_id_fifty.csv \
#  /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/CM/gffcompare_only_r/CM_mando_gene_id_novel_isoforms.txt \
#  --assembly_type filtered

#  python mandalorion/count_nr_novel_nmd_transcripts.py \
#  /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/HUVEC/fifty_nt_result/HUVEC_mando_gene_id_fifty.csv \
#  /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/HUVEC/gffcompare_only_r/HUVEC_mando_gene_id_novel_isoforms.txt \
#  --assembly_type filtered



#################################################################################
# ------------------ COMPARE TO ENSEMBL FULL  ASSEMBLY       ------------------ #
#################################################################################


# if [[ ! -d "$out_path/HUVEC/compare_Ens_full_ref" ]]; then
#     mkdir $out_path/HUVEC/compare_Ens_full_ref
# fi

# if [[ ! -d "$out_path/CM/compare_Ens_full_ref" ]]; then
#     mkdir $out_path/CM/compare_Ens_full_ref
# fi

# if [[ ! -d "$out_path/HUVEC/compare_Ens_full_ref/gffcompare_only_r" ]]; then
#     mkdir $out_path/HUVEC/compare_Ens_full_ref/gffcompare_only_r
# fi

# if [[ ! -d "$out_path/CM/compare_Ens_full_ref/gffcompare_only_r" ]]; then
#     mkdir $out_path/CM/compare_Ens_full_ref/gffcompare_only_r
# fi



# gffcompare -o $out_path/HUVEC/compare_Ens_full_ref/gffcompare_only_r/HUVEC_compare_full_GTF\
#  -r $ensembl_full_gtf\
#   $out_path/HUVEC/HUVEC_mando_gene_id.gtf

# mv $out_path/HUVEC/HUVEC_compare_full_GTF* $out_path/HUVEC/compare_Ens_full_ref/gffcompare_only_r/


# gffcompare -o $out_path/CM/compare_Ens_full_ref/gffcompare_only_r/CM_compare_full_GTF\
#  -r $ensembl_full_gtf\
#   $out_path/CM/CM_mando_gene_id.gtf

# mv $out_path/CM/CM_compare_full_GTF* $out_path/CM/compare_Ens_full_ref/gffcompare_only_r/

# # which isoforms have non ejcs?
# python mandalorion/get_equal_ejc_isoforms.py \
#  $out_path/CM/compare_Ens_full_ref/gffcompare_only_r/CM_compare_full_GTF.CM_mando_gene_id.gtf.tmap

# python mandalorion/get_equal_ejc_isoforms.py \
#  $out_path/HUVEC/compare_Ens_full_ref/gffcompare_only_r/HUVEC_compare_full_GTF.HUVEC_mando_gene_id.gtf.tmap








# which isoforms are novel nmd transcripts?
#  python mandalorion/count_nr_novel_nmd_transcripts.py \
#  /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/CM/fifty_nt_result/CM_mando_gene_id_fifty.csv \
#  $out_path/CM/compare_Ens_full_ref/gffcompare_only_r/CM_mando_gene_id_novel_isoforms.txt \
#  --assembly_type full

#  python mandalorion/count_nr_novel_nmd_transcripts.py \
#  /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/HUVEC/fifty_nt_result/HUVEC_mando_gene_id_fifty.csv \
#  $out_path/HUVEC/compare_Ens_full_ref/gffcompare_only_r/HUVEC_mando_gene_id_novel_isoforms.txt \
#  --assembly_type full


# not sure...maybe keep this step on local to avoid github confusion...
