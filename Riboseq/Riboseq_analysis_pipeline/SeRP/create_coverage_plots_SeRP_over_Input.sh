#!/bin/bash

UMI_dedup_outdir_transcriptomic=$1

# # first merge replicates
# samtools merge -o $UMI_dedup_outdir_transcriptomic/filtered/In_Puro_dedup_q10_filtered_merged.bam\
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_07_In_Puro_1_dedup_q10_filtered.bam\
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_08_In_Puro_3_dedup_q10_filtered.bam \
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_09_In_Puro_4_dedup_q10_filtered.bam 


# samtools merge -o $UMI_dedup_outdir_transcriptomic/filtered/In_CHX_dedup_q10_filtered_merged.bam\
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_10_In_CHX_1_dedup_q10_filtered.bam \
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_11_In_CHX_2_dedup_q10_filtered.bam \
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_12_In_CHX_4_dedup_q10_filtered.bam  

#  samtools merge -o $UMI_dedup_outdir_transcriptomic/filtered/IP_Puro_dedup_q10_filtered_merged.bam\
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_13_IP_Puro_1_dedup_q10_filtered.bam\
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_14_IP_Puro_3_dedup_q10_filtered.bam \
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_15_IP_Puro_4_dedup_q10_filtered.bam

#  samtools merge -o $UMI_dedup_outdir_transcriptomic/filtered/IP_CHX_dedup_q10_filtered_merged.bam\
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_16_IP_CHX_1_dedup_q10_filtered.bam\
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_17_IP_CHX_2_dedup_q10_filtered.bam \
#  $UMI_dedup_outdir_transcriptomic/filtered/uf_muellermcnicoll_2025_04_18_IP_CHX_4_dedup_q10_filtered.bam 

# for bam in "$UMI_dedup_outdir_transcriptomic"/filtered/*merged.bam
# do
#     sorted_bam=$UMI_dedup_outdir_transcriptomic/filtered/$(basename $bam .bam)_sorted.bam
#     samtools sort $bam -o $sorted_bam
#     samtools index $sorted_bam
# done


# bamCompare -b1 $UMI_dedup_outdir_transcriptomic/filtered/IP_CHX_dedup_q10_filtered_merged_sorted.bam\
#  -b2 $UMI_dedup_outdir_transcriptomic/filtered/In_CHX_dedup_q10_filtered_merged_sorted.bam\
#  -o $UMI_dedup_outdir_transcriptomic/filtered/CHX_SeRP_over_Input.bw\
#  -of bigwig\
#  -p 16

# bamCompare -b1 $UMI_dedup_outdir_transcriptomic/filtered/IP_Puro_dedup_q10_filtered_merged_sorted.bam\
#  -b2 $UMI_dedup_outdir_transcriptomic/filtered/In_Puro_dedup_q10_filtered_merged_sorted.bam\
#  -o $UMI_dedup_outdir_transcriptomic/filtered/Puro_SeRP_over_Input.bw\
#  -of bigwig\
#  -p 16


SRSF_MANE_transcripts=( 
"ENST00000258962.5" 
"ENST00000359995.10"
"ENST00000373715.11" 
"ENST00000373795.7" 
"ENST00000557154.6" 
"ENST00000244020.5" 
"ENST00000313117.11"
"ENST00000587424.3"
"ENST00000229390.8"
"ENST00000492112.3"
"ENST00000370949.2"
"ENST00000452027.3" )

# need to get the lengths of the transcripts
# /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai
# has the transcript name as the first column and the length in the second column

for MANE_trans in "${SRSF_MANE_transcripts[@]}"
do
    trans_length=$(grep $MANE_trans  /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai | cut -f2)
    pyGenomeTracks --tracks $UMI_dedup_outdir_transcriptomic/filtered/tracks_CHX_over_SeRP.ini\
    --region $MANE_trans:0-$trans_length\
    --outFileName $UMI_dedup_outdir_transcriptomic/filtered/${MANE_trans}_CHX_over_SeRP.pdf

    pyGenomeTracks --tracks $UMI_dedup_outdir_transcriptomic/filtered/tracks_Puro_over_SeRP.ini\
    --region $MANE_trans:0-$trans_length\
    --outFileName $UMI_dedup_outdir_transcriptomic/filtered/${MANE_trans}_Puro_over_SeRP.pdf

done


