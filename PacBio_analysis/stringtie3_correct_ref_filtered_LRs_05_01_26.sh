#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio


# Define the usage function
usage() {
  echo "Usage: $0 -r <reference_gtf> -g <genome_fasta> -o <out_path> -c <cell_type> -e <ensembl_full_gtf> -n <consensus_reads_fofn> -s <script_dir> -p <short_read_dir>[-h for help]"
}

# Process options with silent error mode
while getopts "b:r:e:f:o:c:n:l:s:p:h" opt; do
  case $opt in
    b)
      bam_dir="$OPTARG"
      ;;
    r)
      reference_gtf="$OPTARG"
      ;;
    e)
     ensembl_full_gtf="$OPTARG"
      ;;
    f)
      genome_fasta="$OPTARG"
      ;;
    o)
      out_path="$OPTARG"
      ;;
    c)
      cell_type="$OPTARG"
      ;;
    n)
      consensus_reads_fofn=="$OPTARG"
      ;;
    l)
      long_read_dir="$OPTARG"
      ;;
    s)
      script_dir="$OPTARG"
      ;;
    p)
      short_read_dir="$OPTARG"
      ;;
    h)
      usage
      exit 0
      ;;
    :)
      echo "Error: Option -$OPTARG requires an argument."
      usage
      exit 1
      ;;
    \?)
      echo "Error: Invalid option -$OPTARG"
      usage
      exit 1
      ;;
  esac
done




if [[ ! -d "$out_path" ]]; then
    mkdir $out_path
fi

if [[ ! -d "$out_path/${cell_type}" ]]; then
    mkdir $out_path/${cell_type}
fi



if [[ ! -d "$bam_dir/merged" ]]; then
    mkdir $bam_dir/merged
fi




#################################################################################
# ------------------ ALIGN LRs IF NECESSARY                      -------------- #
#################################################################################



if [[ ! -d "${bam_dir}" ]]; then
    bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/map_conditions/genome_mapping_cell_type.sh \
    -o "/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment" \
    -f "${genome_fasta}" \
    -i "${long_read_dir}" \
    -c "$cell_type"
fi



#################################################################################
# ------------------ MERGE BAM FILES                         ------------------ #
#################################################################################

# if [ ! -e "$bam_dir/merged/cm_merged.bam" ]; then
#     samtools merge -@ 32 -o $bam_dir/merged/cm_merged.bam \
#     $bam_dir/CM_5NMD_pbmm2_aligned_genome_sorted.bam \
#     $bam_dir/CM_DHYPO_pbmm2_aligned_genome_sorted.bam \
#     $bam_dir/CM_DNOR_pbmm2_aligned_genome_sorted.bam

#     samtools sort -o $bam_dir/merged/cm_merged_sorted.bam $bam_dir/merged/cm_merged.bam

#     samtools index $bam_dir/merged/cm_merged_sorted.bam
# fi

if [ ! -e "$bam_dir/merged/${cell_type}_merged.bam" ]; then
    samtools merge -@ 32 -o $bam_dir/merged/${cell_type}_merged.bam \
    $bam_dir/*filtered.bam

    samtools sort -o $bam_dir/merged/${cell_type}_merged_sorted.bam $bam_dir/merged/${cell_type}_merged.bam

    samtools index $bam_dir/merged/${cell_type}_merged_sorted.bam
fi

#################################################################################
# ------------------ RUN Stringtie   TO CREATE ASSEMBLY      ------------------ #
#################################################################################
if [ ! -e "$out_path/${cell_type}/${cell_type}_strigntie3_assembly.gtf" ]; then
    /home/ckalk/tools/stringtie-3.0.1.Linux_x86_64/stringtie\
    -o $out_path/${cell_type}/${cell_type}_strigntie3_assembly.gtf \
    -L -G $reference_gtf \
    $bam_dir/merged/${cell_type}_merged_sorted.bam
fi


if [[ ! -d "$out_path/${cell_type}/gffcompare" ]]; then
    mkdir $out_path/${cell_type}/gffcompare

    conda activate pacbio
     awk -F " " '$7!="."'  $out_path/${cell_type}/${cell_type}_strigntie3_assembly.gtf\
    >  $out_path/${cell_type}/${cell_type}_strigntie3_assembly_filtered.gtf

    gffcompare -o $out_path/${cell_type}/gffcompare \
    -r $reference_gtf \
    $out_path/${cell_type}/${cell_type}_strigntie3_assembly_filtered.gtf

    conda activate pygtftk
    python /home/ckalk/scripts/SplitORFs/PacBio_analysis/stringtie3/renaming_scripts/rename_STRG_only_genes.py \
    $out_path/${cell_type}/${cell_type}_strigntie3_assembly_filtered.gtf\
    $out_path/${cell_type}/gffcompare.${cell_type}_strigntie3_assembly_filtered.gtf.tmap\
    $out_path/${cell_type}/${cell_type}_strigntie3_assembly_renamed_filtered.gtf
fi



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


if [ ! -d "${stringtie3_dir_raw}"/kallisto/${cell_type}_quant ]; then
    mkdir "${stringtie3_dir_raw}"/kallisto/${cell_type}_quant

    bash ${script_dir}/kallisto/kallisto_index.sh \
    $out_path/${cell_type}/${cell_type}_strigntie3_assembly_filtered.gtf \
    ${genome_fasta} \
    "${stringtie3_dir_raw}"/kallisto/${cell_type}_strigntie3_assembly_transcriptome.fa \
    "${stringtie3_dir_raw}"/kallisto/index/${cell_type}

    bash ${script_dir}/kallisto/kallisto_quantification.sh \
    "${stringtie3_dir_raw}"/kallisto/index/${cell_type}.idx \
    ${short_read_dir} \
    "${stringtie3_dir_raw}"/kallisto/${cell_type}_quant
fi


if [ ! -d "${sqanti_dir}" ]; then
    mkdir "${sqanti_dir}"
fi


if [ ! -d "${sqanti_dir}"/SQANTI3_QC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_QC
fi


if [ ! -d "${sqanti_dir}"/SQANTI3_QC/${cell_type} ]; then
    mkdir "${sqanti_dir}"/SQANTI3_QC/${cell_type}

    bash ${script_dir}/sqanti3/sqanti3_qc_mando_huvec.sh \
    /home/ckalk/tools/sqanti3 \
    $out_path/${cell_type}/${cell_type}_strigntie3_assembly_filtered.gtf \
    ${reference_gtf} \
    ${genome_fasta} \
    "${sqanti_dir}"/SQANTI3_QC/${cell_type} \
    ${script_dir}/sqanti3/${cell_type}_short_reads.txt \
    "${stringtie3_dir_raw}"/kallisto/${cell_type}_quant
fi


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

