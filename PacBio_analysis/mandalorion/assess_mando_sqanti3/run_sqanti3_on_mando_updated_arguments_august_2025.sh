#!/bin/bash

#----- This script performs the preprocessing steps for the short RNA-seq samples ----- #
# ----- and then uses them for SQANTI3 run of the raw HUVEC MAndo assembly        ----- #
# ----- FASTQC is run to obtain quality control metrics and plots                 ----- #
# ----- kallisto for quantification as required for the SQANTI3 pipeline          ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio



################################################################################
# PATH DEFINTIONS                                                              #
################################################################################
# preprocessing directories
# raw_data_dir="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData"
# merged_data_dir="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/merged"
# raw_data_fastqc_dir="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData/fastqc"
# multiQC_outname="multiqc"

# reference file directories
genome_fasta=$1
reference_gtf=$2
mando_out_dir=$3
gtf_file=${mando_out_dir}"/HUVEC/HUVEC_mando_gene_id.gtf"
gtf_file_cm=${mando_out_dir}"/CM/CM_mando_gene_id.gtf"
decoys="/projects/splitorfs/work/reference_files/decoys.txt"
transcript_fasta=${mando_out_dir}"/HUVEC/HUVEC_mando_gene_id_correct.fasta"
transcript_cm_fasta=${mando_out_dir}"/CM/CM_mando_gene_id_correct.fasta"

huvec_lrs=${mando_out_dir}/HUVEC/HUVEC_fl_counts.tsv
cm_lrs=${mando_out_dir}/CM/CM_fl_counts.tsv

mando_dir_raw=$4
sqanti_dir=${mando_out_dir}"/SQANTI3"
outdir_fastp=$5
outdir_fastp_cm=$7
kallisto_quant_mando_raw=${mando_dir_raw}"/kallisto/quant"
kallisto_quant_mando_raw_cm=${mando_dir_raw}"/kallisto/quant_cm"
kallisto_index_path="${mando_dir_raw}"/kallisto/index/HUVEC
kallisto_index_cm_path="${mando_dir_raw}"/kallisto/index/CM
sqanti_qc_outdir="${sqanti_dir}"/SQANTI3_QC

script_dir=$6

isoseq_reads_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"





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


if [ ! -d "${mando_dir_raw}"/kallisto/quant ]; then
    mkdir "${mando_dir_raw}"/kallisto/quant
fi

if [ ! -d "${mando_dir_raw}"/kallisto/quant_cm ]; then
    mkdir "${mando_dir_raw}"/kallisto/quant_cm
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

if [ ! -d "${sqanti_dir}"/SQANTI3_Filter/HUVEC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Filter/HUVEC
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/HUVEC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/HUVEC
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/HUVEC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/CM
fi


################################################################################
# PREPROCESSING STEPS                                                          #
################################################################################
# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/fastqc_multiqc_RNA_seq.sh ${raw_data_dir} ${raw_data_fastqc_dir} ${multiQC_outname} raw

# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/merge_fastq_files.sh ${raw_data_dir} ${merged_data_dir}

# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/fastqc_multiqc_RNA_seq.sh ${merged_data_dir} ${merged_data_dir}/fastqc ${multiQC_outname} raw

# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/trim_adapters_RNA_seq.sh ${merged_data_dir} ${outidr_fastp} _fastp

# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/fastqc_multiqc_RNA_seq.sh ${outidr_fastp} ${outidr_fastp}/fastqc ${multiQC_outname} fastp


################################################################################
# KALLISTO ON RAW MANDO HUVEC ASSEMBLY                                         #
################################################################################
# bash ${script_dir}/kallisto/kallisto_index.sh \
#   ${gtf_file} \
#   ${genome_fasta} \
#   ${transcript_fasta} \
#   ${kallisto_index_path}

# bash ${script_dir}/kallisto/kallisto_quantification.sh \
#  ${kallisto_index_path}.idx \
#  ${outdir_fastp} \
#  ${kallisto_quant_mando_raw}


bash ${script_dir}/kallisto/kallisto_index.sh \
  ${gtf_file_cm} \
  ${genome_fasta} \
  ${transcript_cm_fasta} \
  ${kallisto_index_cm_path}

bash ${script_dir}/kallisto/kallisto_quantification.sh \
 ${kallisto_index_cm_path}.idx \
 ${outdir_fastp_cm} \
 ${kallisto_quant_mando_raw_cm}


################################################################################
# SQANTI3 QC on HUVEC and CM assembly                                          #
################################################################################
 # get the txt file of R1 space R2
# short_read_file=${script_dir}/"sqanti3/huvec_short_reads.txt"

# > "$short_read_file"  # Clear or create the output file

# for r1 in $outdir_fastp/*_merged_fastp.R1.fastp.fastq.gz; do
#     r2="${r1/R1/R2}"
#     if [ -f "$r2" ]; then
#         echo "$r1 $r2" >> "$short_read_file"
#     else
#         echo "Warning: No matching R2 for $r1" >&2
#     fi
# done

# if [ ! -d "${sqanti_qc_outdir}"/CM ]; then
#     mkdir "${sqanti_qc_outdir}"/CM
# fi

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  ${gtf_file_cm} \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  ${sqanti_qc_outdir}/CM \
#  ${cm_lrs}

# if [ ! -d "${sqanti_qc_outdir}"/HUVEC ]; then
#     mkdir "${sqanti_qc_outdir}"/HUVEC
# fi

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_huvec.sh \
#  /home/ckalk/tools/sqanti3 \
#  ${gtf_file} \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  ${sqanti_qc_outdir}/HUVEC \
#  ${short_read_file} \
#  ${kallisto_quant_mando_raw} \
#  ${huvec_lrs}






################################################################################
# SQANTI3 RULES FILTER on HUVEC assembly                                       #
################################################################################
bash ${script_dir}/sqanti3/sqanti_rules/sqanti3_rules_06_08_25.sh \
 /home/ckalk/tools/sqanti3 \
 ${sqanti_qc_outdir}/HUVEC/isoforms \
 ${sqanti_dir}/SQANTI3_Filter/HUVEC \
 ${script_dir}/sqanti3/sqanti_rules/logic_filter_v5_02_09_25.json 

bash ${script_dir}/sqanti3/sqanti_rules/sqanti3_rules_06_08_25.sh \
 /home/ckalk/tools/sqanti3 \
 ${sqanti_qc_outdir}/CM/isoforms \
 ${sqanti_dir}/SQANTI3_Filter/CM \
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
    mkdir "${reference_dir}"/kallisto/quant
fi


kallisto_quant_mando_raw_reference=${reference_dir}/kallisto


################################################################################
# KALLISTO ON RAW ENSEMBL FILTERED REFERENCE                                   #
################################################################################

# bash ${script_dir}/kallisto/kallisto_index.sh \
#   ${reference_gtf} \
#   ${genome_fasta} \
#   ${kallisto_quant_mando_raw_reference}/Ens_110_filtered_transcriptome.fa \
#   ${kallisto_quant_mando_raw_reference}/index

# bash ${script_dir}/kallisto/kallisto_quantification.sh \
#  ${kallisto_quant_mando_raw_reference}/index.idx \
#  ${outdir_fastp} \
#  ${kallisto_quant_mando_raw_reference}/quant


################################################################################
# GET FL COUNTS FOR REFERENCE                                                  #
################################################################################
if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered
fi

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/HUVEC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/HUVEC
fi
if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/CM ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/CM
fi

# pbmm2 index ${kallisto_quant_mando_raw_reference}/Ens_110_filtered_transcriptome.fa \
#  "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/Ens_index.mmi --preset ISOSEQ

# for bam in "${isoseq_reads_dir}"/HUVEC*bam; do
#     sample=$(basename $bam)
#     sample="${sample%_merged_lima_refined.bam}"
#     pbmm2 align "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/Ens_index.mmi \
#     $bam \
#     "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/HUVEC/${sample}_pbmm2_aligned.bam \
#     --preset ISOSEQ \
#     --log-level INFO

#     bam="${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/HUVEC/${sample}_pbmm2_aligned.bam
#     samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
#     samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
#     samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_idxstats.out
#      samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_flagstat.out
#  done

#  for bam in "${isoseq_reads_dir}"/CM*bam; do
#     sample=$(basename $bam)
#     sample="${sample%_merged_lima_refined.bam}"
#     pbmm2 align "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/Ens_index.mmi \
#     $bam \
#     "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/CM/${sample}_pbmm2_aligned.bam \
#     --preset ISOSEQ \
#     --log-level INFO

#     bam="${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/CM/${sample}_pbmm2_aligned.bam
#     samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
#     samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
#     samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_idxstats.out
#      samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_flagstat.out
#  done

# python ${script_dir}/sqanti3/sqanti_rescue/idxstats_for_fl_counts.py \
#  "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/CM

# python ${script_dir}/sqanti3/sqanti_rescue/idxstats_for_fl_counts.py \
#  "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/HUVEC


# SR_bam_file_reference=${script_dir}/"sqanti3/sqanti_rescue/SR_bam_Ens_110.fofn"

# > "$SR_bam_file_reference"  # Clear or create the output file

# for bam in ${reference_dir}/STAR/*sorted.bam; do
#     echo "$bam" >> "$SR_bam_file_reference"
# done




# if [ ! -d ${sqanti_qc_outdir}/Ens_110_filtered_QC/HUVEC ]; then
#     mkdir ${sqanti_qc_outdir}/Ens_110_filtered_QC/HUVEC
# fi
# if [ ! -d "${sqanti_qc_outdir}/Ens_110_filtered_QC/CM" ]; then
#     mkdir ${sqanti_qc_outdir}/Ens_110_filtered_QC/CM
# fi

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_ref_huvec.sh \
#  /home/ckalk/tools/sqanti3 \
#  ${reference_gtf}\
#  ${reference_gtf} \
#  ${genome_fasta} \
#  ${sqanti_qc_outdir}/Ens_110_filtered_QC/HUVEC \
#  ${kallisto_quant_mando_raw_reference} \
#  "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/HUVEC/HUVEC_idx_fl_counts.txt \
#  $SR_bam_file_reference \
#  ${reference_dir}/STAR

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  ${reference_gtf} \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  ${sqanti_qc_outdir}/Ens_110_filtered_QC/CM \
#  "${sqanti_dir}"/SQANTI3_Rescue/pbmm2_Ens_filtered/CM/CM_idx_fl_counts.txt


################################################################################
# RUN THE RESCUE                                                               #
################################################################################


bash ${script_dir}/sqanti3/sqanti_rescue/sqanti_rescue.sh \
    /home/ckalk/tools/sqanti3 \
    ${sqanti_dir}/SQANTI3_Filter/HUVEC/isoforms_classification_TPM.tx.filtered.gtf \
    /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_QC/Ens_110_filtered_QC/HUVEC/isoforms_corrected.gtf \
    ${genome_fasta} \
    ${sqanti_dir}/SQANTI3_Filter/HUVEC/isoforms_classification_TPM.tx_RulesFilter_result_classification.txt \
    ${script_dir}/sqanti3/sqanti_rules/logic_filter_v5_02_09_25.json \
    HUVEC_rescue_rules_filter \
    ${sqanti_qc_outdir}/Ens_110_filtered_QC/HUVEC/isoforms_classification_TPM.txt \
    ${sqanti_qc_outdir}/HUVEC/isoforms_corrected.fasta \
    "${sqanti_dir}"/SQANTI3_Rescue/HUVEC


bash ${script_dir}/sqanti3/sqanti_rescue/sqanti_rescue.sh \
    /home/ckalk/tools/sqanti3 \
    ${sqanti_dir}/SQANTI3_Filter/CM/isoforms_classification_TPM.tx.filtered.gtf \
    /projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_QC/Ens_110_filtered_QC/CM/isoforms_corrected.gtf \
    ${genome_fasta} \
    ${sqanti_dir}/SQANTI3_Filter/CM/isoforms_classification_TPM.tx_RulesFilter_result_classification.txt \
    ${script_dir}/sqanti3/sqanti_rules/logic_filter_v5_02_09_25.json \
    CM_rescue_rules_filter \
    ${sqanti_qc_outdir}/Ens_110_filtered_QC/CM/isoforms_classification_TPM.txt \
    ${sqanti_qc_outdir}/CM/isoforms_corrected.fasta \
    "${sqanti_dir}"/SQANTI3_Rescue/CM


if [ ! -d ${sqanti_dir}/SQANTI3_Rescue/CM/QC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/CM/QC
fi

bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
 /home/ckalk/tools/sqanti3 \
 "${sqanti_dir}"/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf \
 ${reference_gtf} \
 ${genome_fasta} \
 "${sqanti_dir}"/SQANTI3_Rescue/CM/QC

if [ ! -d "${sqanti_dir}"/SQANTI3_Rescue/HUVEC/QC ]; then
    mkdir "${sqanti_dir}"/SQANTI3_Rescue/HUVEC/QC
fi

bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
 /home/ckalk/tools/sqanti3 \
 "${sqanti_dir}"/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf \
 ${reference_gtf} \
 ${genome_fasta} \
 "${sqanti_dir}"/SQANTI3_Rescue/HUVEC/QC





