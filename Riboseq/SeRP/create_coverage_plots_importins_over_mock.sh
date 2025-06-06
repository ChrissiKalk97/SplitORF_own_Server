#!/bin/bash

bam_dir=$1

# first merge replicates
# samtools merge -b <(ls  $bam_dir/*RR_M_CHX*.bam) -o  $bam_dir/enrichment_plots/RR_M_CHX_merged.bam
# samtools merge -b <(ls  $bam_dir/*RR_In_CHX*.bam) -o  $bam_dir/enrichment_plots/RR_In_CHX_merged.bam
# samtools merge -b <(ls  $bam_dir/*RR_A2_CHX*.bam) -o  $bam_dir/enrichment_plots/RR_A2_CHX_merged.bam
# samtools merge -b <(ls  $bam_dir/*RR_B1_CHX*.bam) -o  $bam_dir/enrichment_plots/RR_B1_CHX_merged.bam

# samtools merge -b <(ls  $bam_dir/*RR_A2_Puro*.bam) -o  $bam_dir/enrichment_plots/RR_A2_Puro_merged.bam
# samtools merge -b <(ls  $bam_dir/*RR_B1_Puro*.bam) -o  $bam_dir/enrichment_plots/RR_B1_Puro_merged.bam
# samtools merge -b <(ls  $bam_dir/*RR_In_Puro*.bam) -o  $bam_dir/enrichment_plots/RR_In_Puro_merged.bam


# for bam in "$bam_dir"/enrichment_plots/*merged.bam
# do
#     sorted_bam=$bam_dir/enrichment_plots/$(basename $bam .bam)_sorted.bam
#     samtools sort $bam -o $sorted_bam
#     samtools index $sorted_bam
# done


# bamCompare -b1 $bam_dir/enrichment_plots/RR_A2_CHX_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_M_CHX_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/CHX_A2_over_Mock.bw\
#  -of bigwig\
#  -p 16

# bamCompare -b1 $bam_dir/enrichment_plots/RR_B1_CHX_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_M_CHX_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/CHX_B1_over_Mock.bw\
#  -of bigwig\
#  -p 16

bamCompare -b1 $bam_dir/enrichment_plots/RR_A2_CHX_merged_sorted.bam\
 -b2 $bam_dir/enrichment_plots/RR_In_CHX_merged_sorted.bam\
 -o $bam_dir/enrichment_plots/CHX_A2_over_Input.bw\
 -of bigwig\
 -p 16

bamCompare -b1 $bam_dir/enrichment_plots/RR_B1_CHX_merged_sorted.bam\
 -b2 $bam_dir/enrichment_plots/RR_In_CHX_merged_sorted.bam\
 -o $bam_dir/enrichment_plots/CHX_B1_over_Input.bw\
 -of bigwig\
 -p 16

# bamCompare -b1 $bam_dir/enrichment_plots/RR_A2_Puro_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_In_Puro_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/Puro_A2_over_Input.bw\
#  -of bigwig\
#  -p 16

# bamCompare -b1 $bam_dir/enrichment_plots/RR_B1_Puro_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_In_Puro_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/Puro_B1_over_Input.bw\
#  -of bigwig\
#  -p 16

python map_gids_to_MANE_tids.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $bam_dir/enrichment_plots/2025-06-04_RIP-Seq_hits.csv  \
 $bam_dir/enrichment_plots/RIP_hits_MANE_tIDs.txt \
 $bam_dir/enrichment_plots/RIP_hits_gids_to_MANE_tIDs.csv

readarray -t SRSF_MANE_transcripts < $bam_dir/enrichment_plots/RIP_hits_MANE_tIDs.txt


# # need to get the lengths of the transcripts
# # /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai
# # has the transcript name as the first column and the length in the second column

for MANE_trans in "${SRSF_MANE_transcripts[@]}"
do
#     trans_length=$(grep $MANE_trans  /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai | cut -f2)
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_over_mock.ini\
#     --region $MANE_trans:0-$trans_length\
#     --outFileName $bam_dir/enrichment_plots/${MANE_trans}_CHX_Impa_over_Mock.pdf

#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_over_mock.ini\
#     --region $MANE_trans:0-$trans_length\
#     --outFileName $bam_dir/enrichment_plots/${MANE_trans}_CHX_Impb_over_Mock.pdf

#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_puro_over_in.ini\
#     --region $MANE_trans:0-$trans_length\
#     --outFileName $bam_dir/enrichment_plots/${MANE_trans}_Puro_Impa_over_Input.pdf

#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_puro_over_in.ini\
#     --region $MANE_trans:0-$trans_length\
#     --outFileName $bam_dir/enrichment_plots/${MANE_trans}_Puro_Impb_over_Input.pdf


pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_in.ini\
    --region $MANE_trans:0-$trans_length\
    --outFileName $bam_dir/enrichment_plots/${MANE_trans}_CHX_Impa_over_Input.pdf

pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_CHX_over_in.ini\
    --region $MANE_trans:0-$trans_length\
    --outFileName $bam_dir/enrichment_plots/${MANE_trans}_CHX_Impb_over_Input.pdf


done


