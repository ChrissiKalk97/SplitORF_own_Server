#!/bin/bash

#----- This script runs tama merge on 2 assemblies to merge them into one GTF file ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio


# Define the usage function
usage() {
  echo "Usage: $0 -r <reference_gtf> -g <genome_fasta> -o <outdir_tama> -s <stringtie_gtf> -m <mando_gtf> -c <cell_type> -t <tama_tool_path> [-h for help]"
}

# Process options with silent error mode
while getopts "r:f:o:s:m:p:c:t:h" opt; do
  case $opt in
    r)
      reference_gtf="$OPTARG"
      ;;
    f)
      genome_fasta="$OPTARG"
      ;;
    o)
      outdir_tama="$OPTARG"
      ;;
    s)
      stringtie_gtf="$OPTARG"
      ;;
    m)
      mando_gtf="$OPTARG"
      ;;
    p)
      outdir_fastp="$OPTARG"
      ;;
    c)
      cell_type="$OPTARG"
      ;;
    t)
      tama_tool_path="$OPTARG"
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




# reference_gtf="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"
# genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

# mando_rescued_cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf"
# stringtie_cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3/CM/CM_strigntie3_assembly_renamed_filtered.gtf"

# mando_rescued_huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf"
# stringtie_huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3/HUVEC/HUVEC_strigntie3_assembly_renamed_filtered.gtf"

# outdir_tama="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama"

# tama_tool_path="/home/ckalk/tools/tama"

#################################################################################
# ------------------ MERGE WITH TAMA                          ----------------- #
#################################################################################

cell_type="${cell_type^^}"
cell_type_small="${cell_type,,}"

outdir_fastp=${outdir_fastp}/${cell_type}_fastp

if [[ ! -d "$outdir_tama" ]]; then
    mkdir $outdir_tama
fi

if [[ ! -d "$outdir_tama/${cell_type}" ]]; then
    mkdir $outdir_tama/${cell_type}

    # Need to change the order of gene id and transcript ID
    # these are gene_id; transcript_id but need to be transcript_id,; gene_id in 9th field of GTF
    mando_dir=$(dirname $mando_gtf)
    mando_base=$(basename $mando_gtf .gtf)
    mando_tama_gtf=${mando_dir}/${mando_base}_for_tama.gtf

    python /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/change_order_gene_id_transcript_id_sqanti_gtf.py \
    $mando_gtf \
    $mando_tama_gtf



    # Need bed12 format for tama, convert with tama scripts Mando and Stringtie3
    conda activate tama

    # convert Mando CM assembly to BED12
    python ${tama_tool_path}/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py \
    $mando_tama_gtf \
    $mando_dir/$mando_base.bed12


    # convert Stringtie3 CM assembly to BED12
    stringtie_dir=$(dirname $stringtie_gtf)
    stringtie_base=$(basename $stringtie_gtf .gtf)
    python ${tama_tool_path}/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py \
    $stringtie_gtf \
    $stringtie_dir/$stringtie_base.bed12

    # add the respective samples in the file_list.txt for tama
    sed -i "s|annotation_capped.bed|$stringtie_dir/$stringtie_base.bed12|" /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/file_list_${cell_type_small}.txt
    sed -i "s|annotation_nocap.bed|$mando_dir/$mando_base.bed12|" /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/file_list_${cell_type_small}.txt

    sed -i 's/ \+/\t/g' /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/file_list_${cell_type_small}.txt
    sed -i '/^$/d' /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/file_list_${cell_type_small}.txt


    # 
    python ${tama_tool_path}/tama_merge.py \
    -f /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/file_list_${cell_type_small}.txt \
    -p $outdir_tama/${cell_type}/${cell_type}_merged_tama \
    -s mandalorion \
    -a 50 \
    -m 5 \
    -z 50



    python ${tama_tool_path}/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
    $outdir_tama/${cell_type}/${cell_type}_merged_tama.bed \
    $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf

fi
  
#################################################################################
# ------------------ Kallisto for SQANTI QC                   ----------------- #
#################################################################################


if [ ! -d "$outdir_tama"/kallisto ]; then
    mkdir "$outdir_tama"/kallisto
fi

if [ ! -d "$outdir_tama"/kallisto/index ]; then
    mkdir "$outdir_tama"/kallisto/index
fi

if [ ! -e "$outdir_tama"/kallisto/index/${cell_type}.idx ]; then
    

    bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/kallisto/kallisto_index.sh \
    $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf \
    ${genome_fasta} \
    "$outdir_tama"/kallisto/${cell_type}_tama_merged_assembly_transcriptome.fa \
    "$outdir_tama"/kallisto/index/${cell_type}
fi


if [ ! -d "$outdir_tama/kallisto/${cell_type}_quant" ]; then
    mkdir "$outdir_tama/kallisto/${cell_type}_quant"


    bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/kallisto/kallisto_quantification.sh \
    "$outdir_tama/kallisto/index/${cell_type}.idx" \
    ${outdir_fastp} \
    "$outdir_tama/kallisto/${cell_type}_quant"
fi



#################################################################################
# ------------------              SQANTI QC                   ----------------- #
#################################################################################
if [ ! -d "${outdir_tama}"/SQANTI3_QC ]; then
    mkdir "${outdir_tama}"/SQANTI3_QC
fi

if [ ! -d "${outdir_tama}"/SQANTI3_QC/${cell_type} ]; then
    mkdir "${outdir_tama}"/SQANTI3_QC/${cell_type}

    conda activate pacbio

    bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/sqanti3/sqanti3_qc_mando_huvec.sh \
    /home/ckalk/tools/sqanti3 \
    $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf \
    ${reference_gtf} \
    ${genome_fasta} \
    "${outdir_tama}"/SQANTI3_QC/${cell_type} \
    /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/sqanti3/${cell_type}_short_reads.txt \
    "$outdir_tama/kallisto/${cell_type}_quant"
    

    python /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/get_gene_id_tama_gtf.py \
    $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf \
    $outdir_tama/SQANTI3_QC/${cell_type}/isoforms_classification.txt  \
    $outdir_tama/${cell_type}/${cell_type}_merged_tama_gene_id.gtf

    python /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/add_source_to_tama_gtf.py \
    $outdir_tama/${cell_type}/${cell_type}_merged_tama_gene_id.gtf \
    $outdir_tama/${cell_type}/${cell_type}_merged_tama_trans_report.txt

fi

