#!/bin/bash

# ----- this script runs SQANTI3 on CM or HUVEC mando assembly                   ----- #
# ----- FASTQC is run to obtain quality control metrics and plots                 ----- #
# ----- kallisto for quantification as required for the SQANTI3 pipeline          ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio



################################################################################
# PATH DEFINTIONS                                                              #
################################################################################
# reference file directories
cell_type=$1
genome_fasta=$2
reference_gtf=$3
mando_out_dir=$4
mando_dir_raw=$5
sqanti_dir=${mando_out_dir}"/SQANTI3"
outdir_fastp=$6
script_dir=$7
isoseq_reads_dir=$8



if [ -n "${9}" ]; then
    assembly="${9}"
    if [[ $assembly == "isoquant" ]]; then
        gtf_file="${mando_out_dir}/${cell_type}/${cell_type}.transcript_models.gtf"
        transcript_fasta="${mando_dir_raw}/kallisto/${cell_type}_isoquant_assembly_transcriptome.fa"
        lr_count_tsv="${mando_out_dir}/${cell_type}/${cell_type}.discovered_transcript_grouped_file_name_counts.tsv"

    elif [[ $assembly == "stringtie" ]]; then
        gtf_file="${mando_out_dir}"/"${cell_type}"/"${cell_type}"_strigntie3_assembly_filtered.gtf
        transcript_fasta="${mando_dir_raw}/kallisto/${cell_type}"_strigntie3_assembly_transcriptome.fa
        lr_count_tsv=${mando_out_dir}/${cell_type}/${cell_type}_stringtie_quant/${cell_type}/${cell_type}.transcript_grouped_file_name_counts.tsv
    
    fi
else
        gtf_file=${mando_out_dir}"/${cell_type}/${cell_type}_mando_gene_id.gtf"
        transcript_fasta=${mando_out_dir}"/${cell_type}/${cell_type}_mando_gene_id_correct.fasta"
        lr_count_tsv=${mando_out_dir}/${cell_type}/${cell_type}_fl_counts.tsv
fi


decoys="/projects/splitorfs/work/reference_files/decoys.txt"
kallisto_quant_mando_raw="${mando_dir_raw}/kallisto/quant_${cell_type}"
kallisto_index_path="${mando_dir_raw}"/kallisto/index/${cell_type}
sqanti_qc_outdir="${sqanti_dir}"/SQANTI3_QC





################################################################################
# DIRECTORY CREATION                                                           #
################################################################################
if [ ! -d "${mando_dir_raw}" ]; then
    mkdir "${mando_dir_raw}"
fi

if [ ! -d "${mando_dir_raw}"/kallisto ]; then
    mkdir "${mando_dir_raw}"/kallisto
fi


if [ ! -d "${sqanti_dir}" ]; then
    mkdir "${sqanti_dir}"
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_QC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_QC
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Filter ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Filter
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue
fi


################################################################################
# KALLISTO ON RAW MANDO HUVEC ASSEMBLY                                         #
################################################################################

if [[ ! -f "${kallisto_index_path}.idx" ]]; then
    mkdir -p "${mando_dir_raw}"/kallisto/index

    bash ${script_dir}/kallisto/kallisto_index.sh \
    ${gtf_file} \
    ${genome_fasta} \
    ${transcript_fasta} \
    ${kallisto_index_path}
fi

if [ ! -d "${kallisto_quant_mando_raw}" ]; then
    mkdir "${kallisto_quant_mando_raw}"
    bash ${script_dir}/kallisto/kallisto_quantification.sh \
    ${kallisto_index_path}.idx \
    ${outdir_fastp} \
    ${kallisto_quant_mando_raw}
fi



################################################################################
# SQANTI3 QC on HUVEC and CM assembly                                          #
################################################################################
# get the txt file of R1 space R2, needed for SQANTI short read evidence
short_read_file=${script_dir}/"sqanti3/${cell_type}_short_reads.txt"

> "$short_read_file"  # Clear or create the output file

for r1 in $outdir_fastp/*_merged_fastp.R1.fastp.fastq.gz; do
    r2="${r1/fastp.R1/fastp.R2}"
    if [ -f "$r2" ]; then
        echo "$r1 $r2" >> "$short_read_file"
    else
        echo "Warning: No matching R2 for $r1" >&2
    fi
done

if [ ! -d "${sqanti_qc_outdir}"/${cell_type} ]; then
    mkdir "${sqanti_qc_outdir}"/${cell_type}
    bash ${script_dir}/sqanti3/sqanti3_qc_mando_huvec.sh \
    /home/ckalk/tools/sqanti3.6 \
    ${gtf_file} \
    ${reference_gtf} \
    ${genome_fasta} \
    ${sqanti_qc_outdir}/${cell_type} \
    ${short_read_file} \
    ${kallisto_quant_mando_raw} \
    ${lr_count_tsv}
fi



# ################################################################################
# # SQANTI3 RULES FILTER on HUVEC assembly                                       #
# ################################################################################
# accidentially overwrite the old logic filter
if [ ! -d "${sqanti_dir}"/SQANTI3_Filter/${cell_type} ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Filter/${cell_type}
    bash ${script_dir}/sqanti3/sqanti_rules/sqanti3_rules_06_08_25.sh \
    /home/ckalk/tools/sqanti3.6 \
    ${sqanti_qc_outdir}/${cell_type}/isoforms \
    ${sqanti_dir}/SQANTI3_Filter/${cell_type} \
    ${script_dir}/sqanti3/sqanti_rules/logic_filter_v6_01_06_26.json 
fi




################################################################################
# SQANTI3 RESCUE FILTER on HUVEC assembly                                      #
################################################################################
reference_dir="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Ens_110_filtered"

if [ ! -d "${reference_dir}" ]; then
    mkdir "${reference_dir}"
fi

if [ ! -d "${reference_dir}"/kallisto ]; then
    mkdir "${reference_dir}"/kallisto
fi


kallisto_quant_mando_raw_reference="${reference_dir}"/kallisto


# ###############################################################################
# KALLISTO ON RAW ENSEMBL FILTERED REFERENCE                                   #
# ###############################################################################

if [[ ! -f "${kallisto_quant_mando_raw_reference}"/index.idx ]]; then
    bash ${script_dir}/kallisto/kallisto_index.sh \
    ${reference_gtf} \
    ${genome_fasta} \
    "${kallisto_quant_mando_raw_reference}"/Ens_110_filtered_transcriptome.fa \
    "${kallisto_quant_mando_raw_reference}"/index
fi


if [[ ! -d "${kallisto_quant_mando_raw_reference}"/quant_${cell_type} ]]; then
    mkdir "${kallisto_quant_mando_raw_reference}"/quant_${cell_type}
    bash "${script_dir}"/kallisto/kallisto_quantification.sh \
    "${kallisto_quant_mando_raw_reference}"/index.idx \
    "${outdir_fastp}" \
    "${kallisto_quant_mando_raw_reference}"/quant_${cell_type}
fi



################################################################################
# GET FL COUNTS FOR REFERENCE                                                  #
################################################################################
# if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered ]; then
#     mkdir "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered
#     pbmm2 index ${kallisto_quant_mando_raw_reference}/Ens_110_filtered_transcriptome.fa \
#     "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/Ens_index.mmi --preset ISOSEQ
# fi

# if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type} ]; then
#     mkdir "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type}

#      for bam in "${isoseq_reads_dir}"/${cell_type}*bam; do
#         sample=$(basename $bam)
#         sample="${sample%_merged_lima_refined.bam}"
#         pbmm2 align "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/Ens_index.mmi \
#         $bam \
#         "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type}/${sample}_pbmm2_aligned.bam \
#         --preset ISOSEQ \
#         --log-level INFO

#         bam="${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type}/${sample}_pbmm2_aligned.bam
#         samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
#         samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
#         samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#         $(dirname $bam)/$(basename $bam .bam)_idxstats.out
#         samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#         $(dirname $bam)/$(basename $bam .bam)_flagstat.out
#     done

#     python ${script_dir}/sqanti3/sqanti_rescue/idxstats_for_fl_counts.py \
#     "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type}
# fi

# idea: more consistent to quantify with isoquant
shopt -s nullglob
isoseq_reads_dir="/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment/${cell_type}/minimap2_align"
ref_qc_out_path="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3_June_2026_minimap2/SQANTI3/SQANTI3_QC/Ens_110_filtered_QC/${cell_type}"
bams=("${isoseq_reads_dir}"/*sorted.bam)
echo "${bams[@]}"

# order important for the Stringtie run!
mkdir -p ${sqanti_qc_outdir}/Ens_110_filtered_QC/
mkdir -p ${sqanti_qc_outdir}/Ens_110_filtered_QC/${cell_type}
mkdir -p "${ref_qc_out_path}"
mkdir -p "${ref_qc_out_path}"/quant_"${cell_type}"

echo "${ref_qc_out_path}"/quant_"${cell_type}"/"${cell_type}"/"${cell_type}".transcript_grouped_file_name_counts.tsv

if [[ ! -e "${ref_qc_out_path}"/quant_"${cell_type}"/"${cell_type}"/"${cell_type}".transcript_grouped_file_name_counts.tsv ]]; then
  conda activate isoquant
  isoquant \
      --reference "$genome_fasta" \
      --genedb "${reference_gtf}" \
      --no_model_construction \
      --data_type pacbio_ccs \
      --polya_trimmed stranded \
      --bam  "${bams[@]}" \
      --output "${ref_qc_out_path}"/quant_"${cell_type}" \
      --prefix "${cell_type}"
fi

echo "${kallisto_quant_mando_raw_reference}"

echo "${ref_qc_out_path}"/quant_"${cell_type}"/"${cell_type}"/"${cell_type}".transcript_grouped_file_name_counts.tsv
if [ ! -e "${ref_qc_out_path}"/isoforms_classification.txt ]; then
    bash ${script_dir}/sqanti3/sqanti3_qc_mando_huvec.sh \
    /home/ckalk/tools/sqanti3.6 \
    ${reference_gtf} \
    ${reference_gtf} \
    ${genome_fasta} \
    "${ref_qc_out_path}"\
    ${short_read_file} \
    ${kallisto_quant_mando_raw_reference}/quant_${cell_type} \
    "${ref_qc_out_path}"/quant_"${cell_type}"/"${cell_type}"/"${cell_type}".transcript_grouped_file_name_counts.tsv
fi





################################################################################
# RUN THE RESCUE                                                               #
################################################################################

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/${cell_type} ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/${cell_type}
    bash ${script_dir}/sqanti3/sqanti_rescue/sqanti_rescue.sh \
        /home/ckalk/tools/sqanti3.6 \
        ${sqanti_dir}/SQANTI3_Filter/${cell_type}/isoforms_classification_TPM.tx.filtered.gtf \
        ${ref_qc_out_path}/isoforms_corrected.gtf \
        ${genome_fasta} \
        ${sqanti_dir}/SQANTI3_Filter/${cell_type}/isoforms_classification_TPM.tx_RulesFilter_classification.txt \
        ${script_dir}/sqanti3/sqanti_rules/logic_filter_v6_01_06_26.json \
        ${cell_type}_rescue_rules_filter \
        ${ref_qc_out_path}/isoforms_classification_TPM.txt \
        ${sqanti_qc_outdir}/${cell_type}/isoforms_corrected.fasta \
        "${sqanti_dir}"/SQANTI3_Rescue/${cell_type}
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/${cell_type}/QC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/${cell_type}/QC
    bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
    /home/ckalk/tools/sqanti3.6 \
    "${sqanti_dir}"/SQANTI3_Rescue/${cell_type}/${cell_type}_rescue_rules_filter_rescued.gtf \
    ${reference_gtf} \
    ${genome_fasta} \
    "${sqanti_dir}"/SQANTI3_Rescue/${cell_type}/QC
fi







