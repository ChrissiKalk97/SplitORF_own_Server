#!/bin/bash

#----- This script runs tama merge on bed files that come from the split TAMA collapse runs  ----- #
#----- and merges them by sample ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio


# Define the usage function
usage() {
  echo "Usage: $0 -r <reference_gtf> -g <genome_fasta> -o <outdir_tama>  -c <cell_type> -t <tama_tool_path> [-h for help]"
}


# Process options with silent error mode
while getopts "o:c:s:t:h" opt; do
  case $opt in
    o)
      outdir_tama="$OPTARG"
      ;;
    c)
      cell_type="$OPTARG"
      ;;
    s)
      tama_script_path="$OPTARG"
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


#################################################################################
# ------------------ TAMA MERGE                                  -------------- #
#################################################################################
shopt -s nullglob 
# iterate over the different samples

# create list of bed files to merge
merge_file_txt=""${tama_script_path}"/file_list_bed_merge_${cell_type}.txt"

if [ ! -e "${merge_file_txt}" ]; then
    # creates new file, overwrites if already present
    : > "${merge_file_txt}"

    for dir in "${outdir_tama}"/"${cell_type}"/*/; do
        sample_name=$(basename "${dir}")

        for bed in "${dir}"/*merged_tama.bed; do
        
            value_col1="no_cap"
            value_col2="1,1,1"
            value_col3="${sample_name}"

            printf "%s\t%s\t%s\t%s\n" \
                "${bed}" "${value_col1}" "${value_col2}" "${value_col3}" >> "${merge_file_txt}"
        done
    done

    conda activate tama

    python "${tama_tool_path}"/tama_merge.py \
    -f "${merge_file_txt}" \
    -p "${outdir_tama}"/"${cell_type}"/"${cell_type}"_merged_tama \
    -a 50 \
    -m 5 \
    -z 50 \
    -d merge_dup
fi











  
# #################################################################################
# # ------------------ Kallisto for SQANTI QC                   ----------------- #
# #################################################################################


# if [ ! -d "$outdir_tama"/kallisto ]; then
#     mkdir "$outdir_tama"/kallisto
# fi

# if [ ! -d "$outdir_tama"/kallisto/index ]; then
#     mkdir "$outdir_tama"/kallisto/index
# fi

# if [ ! -e "$outdir_tama"/kallisto/index/${cell_type}.idx ]; then
    

#     bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/kallisto/kallisto_index.sh \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf \
#     ${genome_fasta} \
#     "$outdir_tama"/kallisto/${cell_type}_tama_merged_assembly_transcriptome.fa \
#     "$outdir_tama"/kallisto/index/${cell_type}
# fi


# if [ ! -d "$outdir_tama/kallisto/${cell_type}_quant" ]; then
#     mkdir "$outdir_tama/kallisto/${cell_type}_quant"


#     bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/kallisto/kallisto_quantification.sh \
#     "$outdir_tama/kallisto/index/${cell_type}.idx" \
#     ${outdir_fastp} \
#     "$outdir_tama/kallisto/${cell_type}_quant"
# fi



# #################################################################################
# # ------------------              SQANTI QC                   ----------------- #
# #################################################################################
# if [ ! -d "${outdir_tama}"/SQANTI3_QC ]; then
#     mkdir "${outdir_tama}"/SQANTI3_QC
# fi

# if [ ! -d "${outdir_tama}"/SQANTI3_QC/${cell_type} ]; then
#     mkdir "${outdir_tama}"/SQANTI3_QC/${cell_type}

#     conda activate pacbio

#     bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/sqanti3/sqanti3_qc_mando_huvec.sh \
#     /home/ckalk/tools/sqanti3 \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf \
#     ${reference_gtf} \
#     ${genome_fasta} \
#     "${outdir_tama}"/SQANTI3_QC/${cell_type} \
#     /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/sqanti3/${cell_type}_short_reads.txt \
#     "$outdir_tama/kallisto/${cell_type}_quant"
    

#     python /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/get_gene_id_tama_gtf.py \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf \
#     $outdir_tama/SQANTI3_QC/${cell_type}/isoforms_classification.txt  \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama_gene_id.gtf

#     python /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/add_source_to_tama_gtf.py \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama_gene_id.gtf \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama_trans_report.txt

# fi

