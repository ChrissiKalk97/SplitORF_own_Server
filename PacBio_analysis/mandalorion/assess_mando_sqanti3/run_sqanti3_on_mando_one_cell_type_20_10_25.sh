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
gtf_file=${mando_out_dir}"/${cell_type}/${cell_type}_mando_gene_id.gtf"
decoys="/projects/splitorfs/work/reference_files/decoys.txt"
transcript_fasta=${mando_out_dir}"/${cell_type}/${cell_type}_mando_gene_id_correct.fasta"


lr_count_tsv=${mando_out_dir}/${cell_type}/${cell_type}_fl_counts.tsv


mando_dir_raw=$5
sqanti_dir=${mando_out_dir}"/SQANTI3"
outdir_fastp=$6
script_dir=$7

kallisto_quant_mando_raw=${mando_dir_raw}"/kallisto/quant_${cell_type}"
kallisto_index_path="${mando_dir_raw}"/kallisto/index/${cell_type}
sqanti_qc_outdir="${sqanti_dir}"/SQANTI3_QC



isoseq_reads_dir=$8





################################################################################
# DIRECTORY CREATION                                                           #
################################################################################
if [ ! -d "${mando_dir_raw}" ]; then
    mkdir "${mando_dir_raw}"
fi

if [ ! -d "${mando_dir_raw}"/kallisto ]; then
    mkdir "${mando_dir_raw}"/kallisto
fi

if [ ! -d "${mando_dir_raw}"/kallisto/index ]; then
    mkdir "${mando_dir_raw}"/kallisto/index
fi


if [ ! -d "${mando_dir_raw}"/kallisto/quant_${cell_type} ]; then
    mkdir "${mando_dir_raw}"/kallisto/quant_${cell_type}
fi


if [ ! -d "${sqanti_dir}" ]; then
    mkdir "${sqanti_dir}"
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_QC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_QC
fi

if [ ! -d "${sqanti_qc_outdir}"/${cell_type} ]; then
    mkdir "${sqanti_qc_outdir}"/${cell_type}
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Filter ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Filter
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Filter/${cell_type} ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Filter/${cell_type}
fi


if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/${cell_type} ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/${cell_type}
fi


################################################################################
# KALLISTO ON RAW MANDO HUVEC ASSEMBLY                                         #
################################################################################
bash ${script_dir}/kallisto/kallisto_index.sh \
  ${gtf_file} \
  ${genome_fasta} \
  ${transcript_fasta} \
  ${kallisto_index_path}

bash ${script_dir}/kallisto/kallisto_quantification.sh \
 ${kallisto_index_path}.idx \
 ${outdir_fastp} \
 ${kallisto_quant_mando_raw}


################################################################################
# SQANTI3 QC on HUVEC and CM assembly                                          #
################################################################################
# get the txt file of R1 space R2, needed for SQANTI short read evidence
short_read_file=${script_dir}/"sqanti3/${cell_type}_short_reads.txt"

> "$short_read_file"  # Clear or create the output file

for r1 in $outdir_fastp/*_merged_fastp.R1.fastp.fastq.gz; do
    r2="${r1/R1/R2}"
    if [ -f "$r2" ]; then
        echo "$r1 $r2" >> "$short_read_file"
    else
        echo "Warning: No matching R2 for $r1" >&2
    fi
done


bash ${script_dir}/sqanti3/sqanti3_qc_mando_huvec.sh \
 /home/ckalk/tools/sqanti3 \
 ${gtf_file} \
 ${reference_gtf} \
 ${genome_fasta} \
 ${sqanti_qc_outdir}/${cell_type} \
 ${short_read_file} \
 ${kallisto_quant_mando_raw} \
 ${lr_count_tsv}


################################################################################
# SQANTI3 RULES FILTER on HUVEC assembly                                       #
################################################################################
bash ${script_dir}/sqanti3/sqanti_rules/sqanti3_rules_06_08_25.sh \
 /home/ckalk/tools/sqanti3 \
 ${sqanti_qc_outdir}/${cell_type}/isoforms \
 ${sqanti_dir}/SQANTI3_Filter/${cell_type} \
 ${script_dir}/sqanti3/sqanti_rules/logic_filter_v5_02_09_25.json 



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

if [ ! -d "${reference_dir}"/kallisto/index ]; then
    mkdir "${reference_dir}"/kallisto/index
fi


if [ ! -d "${reference_dir}"/kallisto/quant ]; then
    mkdir "${reference_dir}"/kallisto/quant_${cell_type}
fi


kallisto_quant_mando_raw_reference=${reference_dir}/kallisto


################################################################################
# KALLISTO ON RAW ENSEMBL FILTERED REFERENCE                                   #
################################################################################

bash ${script_dir}/kallisto/kallisto_index.sh \
  ${reference_gtf} \
  ${genome_fasta} \
  ${kallisto_quant_mando_raw_reference}/Ens_110_filtered_transcriptome.fa \
  ${kallisto_quant_mando_raw_reference}/index

bash ${script_dir}/kallisto/kallisto_quantification.sh \
 ${kallisto_quant_mando_raw_reference}/index.idx \
 ${outdir_fastp} \
 ${kallisto_quant_mando_raw_reference}/quant_${cell_type}


################################################################################
# GET FL COUNTS FOR REFERENCE                                                  #
################################################################################
if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type} ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type}
fi

pbmm2 index ${kallisto_quant_mando_raw_reference}/Ens_110_filtered_transcriptome.fa \
 "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/Ens_index.mmi --preset ISOSEQ


 for bam in "${isoseq_reads_dir}"/${cell_type}*bam; do
    sample=$(basename $bam)
    sample="${sample%_merged_lima_refined.bam}"
    pbmm2 align "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/Ens_index.mmi \
    $bam \
    "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type}/${sample}_pbmm2_aligned.bam \
    --preset ISOSEQ \
    --log-level INFO

    bam="${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type}/${sample}_pbmm2_aligned.bam
    samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
    samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
    samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
     $(dirname $bam)/$(basename $bam .bam)_idxstats.out
     samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
     $(dirname $bam)/$(basename $bam .bam)_flagstat.out
 done

python ${script_dir}/sqanti3/sqanti_rescue/idxstats_for_fl_counts.py \
 "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type}



if [ ! -d "${sqanti_qc_outdir}/Ens_110_filtered_QC/${cell_type}" ]; then
    mkdir ${sqanti_qc_outdir}/Ens_110_filtered_QC/${cell_type}
fi


bash ${script_dir}/sqanti3/sqanti3_qc_mando_huvec.sh \
 /home/ckalk/tools/sqanti3 \
 ${reference_gtf} \
 ${reference_gtf} \
 ${genome_fasta} \
 ${sqanti_qc_outdir}/Ens_110_filtered_QC/${cell_type} \
 ${short_read_file} \
 ${kallisto_quant_mando_raw_reference} \
 "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/${cell_type}/${cell_type}_idx_fl_counts.txt


################################################################################
# RUN THE RESCUE                                                               #
################################################################################


# bash ${script_dir}/sqanti3/sqanti_rescue/sqanti_rescue.sh \
#     /home/ckalk/tools/sqanti3 \
#     ${sqanti_dir}/SQANTI3_Filter/HUVEC/isoforms_classification_TPM.tx.filtered.gtf \
#     /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_QC/Ens_110_filtered_QC/HUVEC/isoforms_corrected.gtf \
#     ${genome_fasta} \
#     ${sqanti_dir}/SQANTI3_Filter/HUVEC/isoforms_classification_TPM.tx_RulesFilter_result_classification.txt \
#     ${script_dir}/sqanti3/sqanti_rules/logic_filter_v5_02_09_25.json \
#     HUVEC_rescue_rules_filter \
#     ${sqanti_qc_outdir}/Ens_110_filtered_QC/HUVEC/isoforms_classification_TPM.txt \
#     ${sqanti_qc_outdir}/HUVEC/isoforms_corrected.fasta \
#     "${sqanti_dir}"/SQANTI3_Rescue/HUVEC


# bash ${script_dir}/sqanti3/sqanti_rescue/sqanti_rescue.sh \
#     /home/ckalk/tools/sqanti3 \
#     ${sqanti_dir}/SQANTI3_Filter/CM/isoforms_classification_TPM.tx.filtered.gtf \
#     /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_QC/Ens_110_filtered_QC/CM/isoforms_corrected.gtf \
#     ${genome_fasta} \
#     ${sqanti_dir}/SQANTI3_Filter/CM/isoforms_classification_TPM.tx_RulesFilter_result_classification.txt \
#     ${script_dir}/sqanti3/sqanti_rules/logic_filter_v5_02_09_25.json \
#     CM_rescue_rules_filter \
#     ${sqanti_qc_outdir}/Ens_110_filtered_QC/CM/isoforms_classification_TPM.txt \
#     ${sqanti_qc_outdir}/CM/isoforms_corrected.fasta \
#     "${sqanti_dir}"/SQANTI3_Rescue/CM


# if [ ! -d ${sqanti_dir}/SQANTI3_Rescue/CM/QC ]; then
#     mkdir "${sqanti_dir}"/SQANTI3_Rescue/CM/QC
# fi

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  "${sqanti_dir}"/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  "${sqanti_dir}"/SQANTI3_Rescue/CM/QC

# if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/HUVEC/QC ]; then
#     mkdir "${sqanti_dir}"/SQANTI3_Rescue/HUVEC/QC
# fi

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  "${sqanti_dir}"/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  "${sqanti_dir}"/SQANTI3_Rescue/HUVEC/QC





