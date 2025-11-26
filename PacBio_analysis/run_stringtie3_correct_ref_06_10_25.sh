#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

# ensembl_filtered_gtf=$1
# genome_fasta=$2
# consensus_reads_fofn=$3


ensembl_filtered_gtf="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
consensus_reads_fofn_HUVEC="pacbio_consensus_HUVEC.fofn"
consensus_reads_fofn_CM="./pacbio_consensus_CM.fofn"
out_path="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3"
# Would it have been better to align with minimap and use the gtf as a bed file for the splice junctions?
# this adds a bonus to known splice jcts
# I have to say, since Stringite3 anyway uses the reference and the assembly resembles the reference so closely
# I do not think that it matters a lot...
bam_dir_huvec="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/HUVEC/pbmm2_align/genome"
bam_dir_cm="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/CM/pbmm2_align/genome"
script_dir="/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3"
outdir_fastp="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/fastp"

if [[ ! -d "$out_path" ]]; then
    mkdir $out_path
fi

if [[ ! -d "$out_path/HUVEC" ]]; then
    mkdir $out_path/HUVEC
fi

if [[ ! -d "$out_path/CM" ]]; then
    mkdir $out_path/CM
fi

if [[ ! -d "$bam_dir_huvec/merged" ]]; then
    mkdir $bam_dir_huvec/merged
fi

if [[ ! -d "$bam_dir_cm/merged" ]]; then
    mkdir $bam_dir_cm/merged
fi



#################################################################################
# ------------------ MERGE BAM FILES                         ------------------ #
#################################################################################
# samtools merge -@ 32 -o $bam_dir_cm/merged/cm_merged.bam \
#  $bam_dir_cm/CM_5NMD_pbmm2_aligned_genome_sorted.bam \
#  $bam_dir_cm/CM_DHYPO_pbmm2_aligned_genome_sorted.bam \
#  $bam_dir_cm/CM_DNOR_pbmm2_aligned_genome_sorted.bam

# samtools merge -@ 32 -o $bam_dir_huvec/merged/huvec_merged.bam \
#  $bam_dir_huvec/HUVEC_50NMD_pbmm2_aligned_genome_sorted.bam \
#  $bam_dir_huvec/HUVEC_5NMD_pbmm2_aligned_genome_sorted.bam \
#  $bam_dir_huvec/HUVEC_DHYPO_pbmm2_aligned_genome_sorted.bam \
#  $bam_dir_huvec/HUVEC_DMSO_pbmm2_aligned_genome_sorted.bam \
#  $bam_dir_huvec/HUVEC_DNOR_pbmm2_aligned_genome_sorted.bam


# samtools sort -o $bam_dir_cm/merged/cm_merged_sorted.bam $bam_dir_cm/merged/cm_merged.bam
# samtools sort -o $bam_dir_huvec/merged/huvec_merged_sorted.bam $bam_dir_huvec/merged/huvec_merged.bam

# samtools index $bam_dir_cm/merged/cm_merged_sorted.bam
# samtools index $bam_dir_huvec/merged/huvec_merged_sorted.bam

#################################################################################
# ------------------ RUN Stringtie   TO CREATE ASSEMBLY      ------------------ #
#################################################################################
# /home/ckalk/tools/stringtie-3.0.1.Linux_x86_64/stringtie\
#  -o $out_path/HUVEC/HUVEC_strigntie3_assembly.gtf \
#  -L -G $ensembl_filtered_gtf \
#  $bam_dir_huvec/merged/huvec_merged_sorted.bam

 /home/ckalk/tools/stringtie-3.0.1.Linux_x86_64/stringtie\
 -o $out_path/CM/CM_strigntie3_assembly.gtf \
 -L -G $ensembl_filtered_gtf \
 $bam_dir_cm/merged/cm_merged_sorted.bam

if [[ ! -d "$out_path/HUVEC/gffcompare" ]]; then
    mkdir $out_path/HUVEC/gffcompare
fi

if [[ ! -d "$out_path/CM/gffcompare" ]]; then
    mkdir $out_path/CM/gffcompare
fi

conda activate pacbio
# Filter out those entries where the strand is not defined
awk -F " " '$7!="."'  $out_path/CM/CM_strigntie3_assembly.gtf\
 >  $out_path/CM/CM_strigntie3_assembly_filtered.gtf

#  awk -F " " '$7!="."'  $out_path/HUVEC/HUVEC_strigntie3_assembly.gtf\
#  >  $out_path/HUVEC/HUVEC_strigntie3_assembly_filtered.gtf


# gffcompare -o $out_path/HUVEC/gffcompare \
#  -r $ensembl_filtered_gtf \
#   $out_path/HUVEC/HUVEC_strigntie3_assembly_filtered.gtf

gffcompare -o $out_path/CM/gffcompare \
 -r $ensembl_filtered_gtf \
  $out_path/CM/CM_strigntie3_assembly_filtered.gtf

# conda activate pygtftk
# python /home/ckalk/scripts/SplitORFs/PacBio_analysis/stringtie3/renaming_scripts/rename_STRG_only_genes.py \
#  $out_path/HUVEC/HUVEC_strigntie3_assembly_filtered.gtf\
#   $out_path/HUVEC/gffcompare.HUVEC_strigntie3_assembly_filtered.gtf.tmap\
#     $out_path/HUVEC/HUVEC_strigntie3_assembly_renamed_filtered.gtf

python /home/ckalk/scripts/SplitORFs/PacBio_analysis/stringtie3/renaming_scripts/rename_STRG_only_genes.py \
 $out_path/CM/CM_strigntie3_assembly_filtered.gtf\
  $out_path/CM/gffcompare.CM_strigntie3_assembly_filtered.gtf.tmap\
    $out_path/CM/CM_strigntie3_assembly_renamed_filtered.gtf






#################################################################################
# ------------------ RUN SQANTI3 for QC assessment           ------------------ #
#################################################################################
stringtie3_dir_raw="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Stringtie3_raw"
sqanti_dir=$out_path/SQANTI3
if [ ! -d "${stringtie3_dir_raw}" ]; then
    mkdir "${stringtie3_dir_raw}"
fi

if [ ! -d "${stringtie3_dir_raw}"/kallisto ]; then
    mkdir "${stringtie3_dir_raw}"/kallisto
fi

if [ ! -d "${stringtie3_dir_raw}"/kallisto/index ]; then
    mkdir "${stringtie3_dir_raw}"/kallisto/index
fi


if [ ! -d "${stringtie3_dir_raw}"/kallisto/quant ]; then
    mkdir "${stringtie3_dir_raw}"/kallisto/quant
fi

# bash ${script_dir}/kallisto/kallisto_index.sh \
#   $out_path/HUVEC/HUVEC_strigntie3_assembly_filtered.gtf \
#   ${genome_fasta} \
#   "${stringtie3_dir_raw}"/kallisto/HUVEC_strigntie3_assembly_transcriptome.fa \
#   "${stringtie3_dir_raw}"/kallisto/index

# bash ${script_dir}/kallisto/kallisto_quantification.sh \
#  "${stringtie3_dir_raw}"/kallisto/index.idx \
#  ${outdir_fastp} \
#  "${stringtie3_dir_raw}"/kallisto/quant







if [ ! -d "${sqanti_dir}" ]; then
    mkdir "${sqanti_dir}"
fi


if [ ! -d "${sqanti_dir}"/SQANTI3_QC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_QC
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Filter ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Filter
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Filter/HUVEC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Filter/HUVEC
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/HUVEC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/HUVEC
fi


if [ ! -d "${sqanti_dir}"/SQANTI3_QC/CM ]; then
    mkdir "${sqanti_dir}"/SQANTI3_QC/CM
fi

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $out_path/CM/CM_strigntie3_assembly_filtered.gtf \
#  ${ensembl_filtered_gtf} \
#  ${genome_fasta} \
#  "${sqanti_dir}"/SQANTI3_QC/CM

if [ ! -d "${sqanti_dir}"/SQANTI3_QC/HUVEC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_QC/HUVEC
fi

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_huvec.sh \
#  /home/ckalk/tools/sqanti3 \
#  $out_path/HUVEC/HUVEC_strigntie3_assembly_filtered.gtf \
#  ${ensembl_filtered_gtf} \
#  ${genome_fasta} \
#  "${sqanti_dir}"/SQANTI3_QC/HUVEC \
#  ${script_dir}/sqanti3/huvec_short_reads.txt \
#  "${stringtie3_dir_raw}"/kallisto/quant


#################################################################################
# ------------------ RUN SPLIT-ORFs PIPELINE                 ------------------ #
#################################################################################
# bash SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
# $out_path/HUVEC/HUVEC_strigntie3_assembly_renamed_filtered.gtf \
# /home/ckalk/tools/SplitORF_pipeline \
# $genome_fasta


# bash SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
# $out_path/CM/CM_strigntie3_assembly_renamed_filtered.gtf \
# /home/ckalk/tools/SplitORF_pipeline \
# $genome_fasta

#################################################################################
# ------------------ RUN 50nt       PIPELINE                 ------------------ #
#################################################################################
# bash SplitORF_scripts/run_fiftynt_on_assembly.sh \
#     $out_path/CM/CM_strigntie3_assembly_renamed_filtered.gtf \
#     /home/ckalk/tools/NMD_fetaure_composition \
#     $genome_fasta \
#     $ensembl_full_gtf \
#     CM_stringtie_50nt.csv


# bash SplitORF_scripts/run_fiftynt_on_assembly.sh \
#     $out_path/HUVEC/HUVEC_strigntie3_assembly_renamed_filtered.gtf \
#     /home/ckalk/tools/NMD_fetaure_composition \
#     $genome_fasta \
#     $ensembl_full_gtf \
#     HUVEC_stringtie_50nt.csv


#################################################################################
# ------------------ COMPARE TO ENSEMBL FULL  ASSEMBLY       ------------------ #
#################################################################################


# if [[ ! -d "$out_path/HUVEC/compare_Ens_full_ref" ]]; then
#     mkdir $out_path/HUVEC/compare_Ens_full_ref
# fi

# if [[ ! -d "$out_path/CM/compare_Ens_full_ref" ]]; then
#     mkdir $out_path/CM/compare_Ens_full_ref
# fi


# gffcompare -o $out_path/HUVEC/compare_Ens_full_ref/HUVEC_compare_full_GTF\
#  -r $ensembl_full_gtf\
#   $out_path/HUVEC/HUVEC_strigntie3_assembly_renamed_filtered.gtf

# mv $out_path/HUVEC/HUVEC_compare_full_GTF* $out_path/HUVEC/compare_Ens_full_ref


# gffcompare -o $out_path/CM/compare_Ens_full_ref/CM_compare_full_GTF\
#  -r $ensembl_full_gtf\
#   $out_path/CM/CM_strigntie3_assembly_renamed_filtered.gtf

# mv $out_path/CM/CM_compare_full_GTF* $out_path/CM/compare_Ens_full_ref

# # # which isoforms have non ejcs?
# python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/get_equal_ejc_isoforms.py \
#  $out_path/CM/compare_Ens_full_ref/CM_compare_full_GTF.CM_strigntie3_assembly_renamed_filtered.gtf.tmap

# python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/get_equal_ejc_isoforms.py \
#  $out_path/HUVEC/compare_Ens_full_ref/HUVEC_compare_full_GTF.HUVEC_strigntie3_assembly_renamed_filtered.gtf.tmap


# # which isoforms are novel nmd transcripts?
#  python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/count_nr_novel_nmd_transcripts.py \
#  /home/ckalk/tools/NMD_fetaure_composition/Output/CM_stringtie_50nt/CM_stringtie_50nt.csv \
#  $out_path/CM/compare_Ens_full_ref/CM_strigntie3_assembly_renamed_filtered_novel_isoforms.txt \
#  --assembly_type full

#  python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/count_nr_novel_nmd_transcripts.py \
#  /home/ckalk/tools/NMD_fetaure_composition/Output/HUVEC_stringtie_50nt/HUVEC_stringtie_50nt.csv \
#  $out_path/HUVEC/compare_Ens_full_ref/HUVEC_strigntie3_assembly_renamed_filtered_novel_isoforms.txt \
#  --assembly_type full

