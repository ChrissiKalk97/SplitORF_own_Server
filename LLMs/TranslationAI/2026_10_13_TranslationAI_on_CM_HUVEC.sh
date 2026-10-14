#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate translationAI

cm_nmd_transcripts="/home/ckalk/tools/SplitORF_pipeline/Input2023/CM_assembly/CM_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID.fa"
cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/CM/CM_LR_SR_support_filtered.gtf"
cm_trans_indir="/projects/splitorfs/work/LLMs/TranslationAI/Input/CM_assembly"

huvec_nmd_transcripts="/home/ckalk/tools/SplitORF_pipeline/Input2023/HUVEC_assembly/HUVEC_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID.fa"
huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/HUVEC/HUVEC_LR_SR_support_filtered.gtf"
huvec_trans_indir="/projects/splitorfs/work/LLMs/TranslationAI/Input/HUVEC_assembly"

mkdir -p "$cm_trans_indir"
mkdir -p "$huvec_trans_indir"


#################################################################################
# ------------------ PREPARE INPUT DATA                      ------------------ #
#################################################################################

echo "Reformat fasta header"

python /home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/reformat_fasta_header_custom_assembly.py\
 "$cm_nmd_transcripts" \
 "$cm_gtf" \
 "$cm_trans_indir"/CM_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID_for_transAI.fa

python /home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/reformat_fasta_header_custom_assembly.py\
 "$huvec_nmd_transcripts" \
 "$huvec_gtf" \
 "$huvec_trans_indir"/HUVEC_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID_for_transAI.fa


head "$huvec_trans_indir"/HUVEC_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID_for_transAI.fa



#################################################################################
# ------------------ RUNNIG TRANSLATION AI                   ------------------ #
#################################################################################

echo "Running TIS preds with TranslationAI"


# Output is written into Input directory
translationai -I "$huvec_trans_indir"/HUVEC_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID_for_transAI.fa -t 0.0000001,0.0000001
translationai -I "$cm_trans_indir"/CM_NMD_transcripts_ORFanage_FIRST_and_fiftyntrule_pipeline_gID_tID_for_transAI.fa -t 0.0000001,0.0000001


# python /home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/analyze_TIS/analyze_TIS_pipeline.py \
#  /projects/splitorfs/work/LLMs/TranslationAI/Output/NMD_trnascripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt \
#  /projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt


