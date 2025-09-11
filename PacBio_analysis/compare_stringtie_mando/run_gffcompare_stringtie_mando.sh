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



script_dir="/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

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
# ------------------ MERGE WITH STRINGTIE MERGE               ----------------- #
#################################################################################
# outdir_stringtie_merge="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/stringtie_merge"

# if [[ ! -d "$outdir_stringtie_merge" ]]; then
#     mkdir $outdir_stringtie_merge
# fi

# if [[ ! -d "$outdir_stringtie_merge/CM" ]]; then
#     mkdir $outdir_stringtie_merge/CM
# fi

# if [[ ! -d "$outdir_stringtie_merge/HUVEC" ]]; then
#     mkdir $outdir_stringtie_merge/HUVEC
# fi


# /home/ckalk/tools/stringtie-3.0.1.Linux_x86_64/stringtie \
#   --merge \
#   -o $outdir_stringtie_merge/CM/CM_mando_stringtie_merged.gtf\
#     $mando_rescued_cm_gtf \
#     $stringtie_cm_gtf

# /home/ckalk/tools/stringtie-3.0.1.Linux_x86_64/stringtie \
#   --merge \
#   -o $outdir_stringtie_merge/HUVEC/HUVEC_mando_stringtie_merged.gtf \
#     $mando_rescued_huvec_gtf \
#     $stringtie_huvec_gtf



# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_stringtie_merge/CM/CM_mando_stringtie_merged.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_stringtie_merge/CM/SQANTI_QC

#  bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_stringtie_merge/HUVEC/HUVEC_mando_stringtie_merged.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_stringtie_merge/HUVEC/SQANTI_QC

#################################################################################
# ------------------ MERGE WITH TAMA                          ----------------- #
#################################################################################

outdir_tama="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama"

# if [[ ! -d "$outdir_tama" ]]; then
#     mkdir $outdir_tama
# fi

# if [[ ! -d "$outdir_tama/CM" ]]; then
#     mkdir $outdir_tama/CM
# fi

# if [[ ! -d "$outdir_tama/HUVEC" ]]; then
#     mkdir $outdir_tama/HUVEC
# fi


# # Need to change the order of gene id and transcript ID
mando_dir_cm=$(dirname $mando_rescued_cm_gtf)
mando_base_cm=$(basename $mando_rescued_cm_gtf .gtf)
mando_tama_cm_gtf=${mando_dir_cm}/${mando_base_cm}_for_tama.gtf
# python change_order_gene_id_transcript_id_sqanti_gtf.py \
#  $mando_rescued_cm_gtf \
#  $mando_tama_cm_gtf



mando_dir_huvec=$(dirname $mando_rescued_huvec_gtf)
mando_base_huvec=$(basename $mando_rescued_huvec_gtf .gtf)
mando_tama_huvec_gtf=${mando_dir_huvec}/${mando_base_huvec}_for_tama.gtf
# python change_order_gene_id_transcript_id_sqanti_gtf.py \
#  $mando_rescued_huvec_gtf \
#  $mando_tama_huvec_gtf



# # Need bed12 format
# conda activate tama

# python /home/ckalk/tools/tama/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py \
#  $mando_tama_cm_gtf \
#  $mando_dir_cm/$mando_base_cm.bed12


# python /home/ckalk/tools/tama/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py \
#  $mando_tama_huvec_gtf \
#  $mando_dir_huvec/$mando_base_huvec.bed12


# stringtie_dir_cm=$(dirname $stringtie_cm_gtf)
# stringtie_base_cm=$(basename $stringtie_cm_gtf .gtf)
# python /home/ckalk/tools/tama/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py \
#  $stringtie_cm_gtf \
#  $stringtie_dir_cm/$stringtie_base_cm.bed12

# stringtie_dir_huvec=$(dirname $stringtie_huvec_gtf)
# stringtie_base_huvec=$(basename $stringtie_huvec_gtf .gtf)
# python /home/ckalk/tools/tama/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py \
#  $stringtie_huvec_gtf \
#  $stringtie_dir_huvec/$stringtie_base_huvec.bed12


# sed -i "s|annotation_capped.bed|$stringtie_dir_cm/$stringtie_base_cm.bed12|" file_list_cm.txt
# sed -i "s|annotation_nocap.bed|$mando_dir_cm/$mando_base_cm.bed12|" file_list_cm.txt

# sed -i 's/ \+/\t/g' file_list_cm.txt
# sed -i '/^$/d' file_list_cm.txt


# sed -i "s|annotation_capped.bed|$stringtie_dir_huvec/$stringtie_base_huvec.bed12|" file_list_huvec.txt
# sed -i "s|annotation_nocap.bed|$mando_dir_huvec/$mando_base_huvec.bed12|" file_list_huvec.txt

# sed -i 's/ \+/\t/g' file_list_huvec.txt
# sed -i '/^$/d' file_list_huvec.txt


# python /home/ckalk/tools/tama/tama_merge.py \
# -f file_list_cm.txt \
# -p $outdir_tama/CM/CM_merged_tama \
# -s mandalorion \
# -a 50 \
# -m 5 \
# -z 50


# python /home/ckalk/tools/tama/tama_merge.py \
# -f file_list_huvec.txt \
# -p $outdir_tama/HUVEC/HUVEC_merged_tama \
# -s mandalorion \
# -a 50 \
# -m 5 \
# -z 50

# python /home/ckalk/tools/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
#  $outdir_tama/CM/CM_merged_tama.bed \
#  $outdir_tama/CM/CM_merged_tama.gtf

# python /home/ckalk/tools/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
#  $outdir_tama/HUVEC/HUVEC_merged_tama.bed \
#  $outdir_tama/HUVEC/HUVEC_merged_tama.gtf


# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_tama/CM/CM_merged_tama.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_tama/CM/SQANTI_QC

#  bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_tama/HUVEC/HUVEC_merged_tama.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_tama/HUVEC/SQANTI_QC

python get_gene_id_tama_gtf.py \
 $outdir_tama/HUVEC/HUVEC_merged_tama.gtf \
 $outdir_tama/HUVEC/SQANTI_QC/isoforms_classification.txt  \
 $outdir_tama/HUVEC/HUVEC_merged_tama_gene_id.gtf

python get_gene_id_tama_gtf.py \
 $outdir_tama/CM/CM_merged_tama.gtf \
 $outdir_tama/CM/SQANTI_QC/isoforms_classification.txt  \
 $outdir_tama/CM/CM_merged_tama_gene_id.gtf