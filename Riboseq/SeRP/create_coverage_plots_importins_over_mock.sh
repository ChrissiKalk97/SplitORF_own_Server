#!/bin/bash

bam_dir=$1

# first merge replicates
# samtools merge -b <(ls $bam_dir/*RR_M_CHX*.bam) $bam_dir/enrichment_plots/RR_M_CHX_merged.bam
# samtools merge -b <(ls $bam_dir/*RR_In_CHX*.bam) $bam_dir/enrichment_plots/RR_In_CHX_merged.bam
# samtools merge -b <(ls $bam_dir/*RR_A2_CHX*.bam) $bam_dir/enrichment_plots/RR_A2_CHX_merged.bam
# samtools merge -b <(ls $bam_dir/*RR_B1_CHX*.bam) $bam_dir/enrichment_plots/RR_B1_CHX_merged.bam

# samtools merge -b <(ls $bam_dir/*RR_A2_Puro*.bam) $bam_dir/enrichment_plots/RR_A2_Puro_merged.bam
# samtools merge -b <(ls $bam_dir/*RR_B1_Puro*.bam) $bam_dir/enrichment_plots/RR_B1_Puro_merged.bam
# samtools merge -b <(ls $bam_dir/*RR_In_Puro*.bam) $bam_dir/enrichment_plots/RR_In_Puro_merged.bam


# for bam in "$bam_dir"/enrichment_plots/*merged.bam
# do
#     sorted_bam=$bam_dir/enrichment_plots/$(basename $bam .bam)_sorted.bam
#     samtools sort $bam -o $sorted_bam
#     samtools index $sorted_bam
# done



# bamCompare -b1 $bam_dir/enrichment_plots/RR_A2_CHX_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_M_CHX_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/CHX_A2_over_Mock_ratio_mov_avg.bw\
#  --operation ratio \
#  --binSize 1 \
#  --smoothLength 21 \
#  -of bigwig\
#  -p 16

# bamCompare -b1 $bam_dir/enrichment_plots/RR_B1_CHX_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_M_CHX_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/CHX_B1_over_Mock_ratio_mov_avg.bw\
#   --operation ratio \
#  --binSize 1 \
#  --smoothLength 21 \
#  -of bigwig\
#  -p 16

# bamCompare -b1 $bam_dir/enrichment_plots/RR_A2_CHX_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_In_CHX_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/CHX_A2_over_Input_ratio_mov_avg.bw\
#   --operation ratio \
#  --binSize 1 \
#  --smoothLength 21 \
#  -of bigwig\
#  -p 16

# bamCompare -b1 $bam_dir/enrichment_plots/RR_B1_CHX_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_In_CHX_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/CHX_B1_over_Input_ratio_mov_avg.bw\
#   --operation ratio \
#  --binSize 1 \
#  --smoothLength 21 \
#  -of bigwig\
#  -p 16

# bamCompare -b1 $bam_dir/enrichment_plots/RR_A2_Puro_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_In_Puro_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/Puro_A2_over_Input_ratio_mov_avg.bw\
#   --operation ratio \
#  --binSize 1 \
#  --smoothLength 21 \
#  -of bigwig\
#  -p 16

# bamCompare -b1 $bam_dir/enrichment_plots/RR_B1_Puro_merged_sorted.bam\
#  -b2 $bam_dir/enrichment_plots/RR_In_Puro_merged_sorted.bam\
#  -o $bam_dir/enrichment_plots/Puro_B1_over_Input_ratio_mov_avg.bw\
#   --operation ratio \
#  --binSize 1 \
#  --smoothLength 21 \
#  -of bigwig\
#  -p 16

# python map_gids_to_MANE_tids.py \
#  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
#  $bam_dir/enrichment_plots/2025-06-04_RIP-Seq_hits.csv  \
#  $bam_dir/enrichment_plots/RIP_hits_MANE_tIDs.txt \
#  $bam_dir/enrichment_plots/RIP_hits_gids_to_MANE_tIDs.csv \
#  multirow



# # need to get the lengths of the transcripts
# # /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai
# # has the transcript name as the first column and the length in the second column

if [ ! -d $bam_dir/enrichment_plots/impA_DEG_plots ]; then
        mkdir $bam_dir/enrichment_plots/impA_DEG_plots
fi

if [ ! -d $bam_dir/enrichment_plots/impB_DEG_plots ]; then
        mkdir $bam_dir/enrichment_plots/impB_DEG_plots
fi

if [ ! -d $bam_dir/enrichment_plots/hot_candiates ]; then
        mkdir $bam_dir/enrichment_plots/hot_candiates
fi




################################################################################
# Imp a enrichment plots                                                       #
################################################################################

# readarray -t importin_enriched_A2 < /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_A2_CHX_enriched_over_Input_and_Mock_MANE_tIDs.txt

# # get the overall max of ImpA_Mock CHX and set to that scale
# max=$(bigWigInfo $bam_dir/enrichment_plots/CHX_A2_over_Mock_ratio_mov_avg.bw | grep "max:" | awk '{print $2}')
# sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini

# # get the overall max of ImpA_Input CHX and set to that scale
# max=$(bigWigInfo $bam_dir/enrichment_plots/CHX_A2_over_Input_ratio_mov_avg.bw | grep "max:" | awk '{print $2}')
# sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini

# # get the overall max of ImpA_Input Puro and set to that scale
# max=$(bigWigInfo $bam_dir/enrichment_plots/Puro_A2_over_Input_ratio_mov_avg.bw | grep "max:" | awk '{print $2}')
# sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini

# for MANE_trans in "${importin_enriched_A2[@]}"
# do
#     trans_length=$(grep $MANE_trans  /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai | cut -f2)

#     # set number of bin to length of the transcript
#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/impA_DEG_plots/${MANE_trans}_CHX_Impa_over_Mock_ratio.pdf

    
#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/impA_DEG_plots/${MANE_trans}_Puro_Impa_over_Input_ratio.pdf

#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini\
#         --region $MANE_trans:0-$trans_length\
#         --height 40\
#         --outFileName $bam_dir/enrichment_plots/impA_DEG_plots/${MANE_trans}_CHX_Impa_over_Input_ratio.pdf

# done


################################################################################
# Imp b enrichment plots                                                       #
################################################################################

# readarray -t importin_enriched_B1 < /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_B1_CHX_enriched_over_Input_and_Mock_MANE_tIDs.txt

# # get the overall max of ImpA_Mock CHX and set to that scale
# max=$(bigWigInfo $bam_dir/enrichment_plots/CHX_B1_over_Mock_ratio_mov_avg.bw | grep "max:" | awk '{print $2}')
# sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini

# # get the overall max of ImpA_Input CHX and set to that scale
# max=$(bigWigInfo $bam_dir/enrichment_plots/CHX_B1_over_Input_ratio_mov_avg.bw | grep "max:" | awk '{print $2}')
# sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini

# # get the overall max of ImpA_Input Puro and set to that scale
# max=$(bigWigInfo $bam_dir/enrichment_plots/Puro_B1_over_Input_ratio_mov_avg.bw | grep "max:" | awk '{print $2}')
# sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini


# for MANE_trans in "${importin_enriched_B1[@]}"
# do
#     trans_length=$(grep $MANE_trans  /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai | cut -f2)
    
#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/impB_DEG_plots/${MANE_trans}_CHX_Impb_over_Mock_ratio.pdf


#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/impB_DEG_plots/${MANE_trans}_Puro_Impb_over_Input_ratio.pdf


#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini\
#         --region $MANE_trans:0-$trans_length\
#         --height 40\
#         --outFileName $bam_dir/enrichment_plots/impB_DEG_plots/${MANE_trans}_CHX_Impb_over_Input_ratio.pdf

# done



################################################################################
# plot hot candidates on the same scale                                        #
################################################################################
hot_candidates=(ENST00000397885.3 ENST00000330560.8 ENST00000428849.7 ENST00000611405.5)

for MANE_trans in "${hot_candidates[@]}"
do
        trans_length=$(grep $MANE_trans  /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai | cut -f2)
        # max_val_A2_CHX_mock=$(bigWigSummary -type=max $bam_dir/enrichment_plots/CHX_A2_over_Mock_ratio_mov_avg.bw $MANE_trans 0 $trans_length 10)
        # max_val_A2_CHX_Input=$(bigWigSummary -type=max $bam_dir/enrichment_plots/CHX_A2_over_Input_ratio_mov_avg.bw $MANE_trans 0 $trans_length 10)
        # max_val_A2_Puro=$(bigWigSummary -type=max $bam_dir/enrichment_plots/Puro_A2_over_Input_ratio_mov_avg.bw $MANE_trans 0 $trans_length 1)
        # max_val_B1_CHX_mock=$(bigWigSummary -type=max $bam_dir/enrichment_plots/CHX_B1_over_Mock_ratio_mov_avg.bw $MANE_trans 0 $trans_length 10)
        # max_val_B1_CHX_Input=$(bigWigSummary -type=max $bam_dir/enrichment_plots/CHX_B1_over_Input_ratio_mov_avg.bw $MANE_trans 0 $trans_length 10)
        # max_val_B1_Puro=$(bigWigSummary -type=max $bam_dir/enrichment_plots/Puro_B1_over_Input_ratio_mov_avg.bw $MANE_trans 0 $trans_length 1)


        # numbers=($max_val_A2_CHX_mock $max_val_A2_CHX_Input $max_val_B1_CHX_mock $max_val_B1_CHX_Input)
        # numbers=($max_val_A2_CHX_mock $max_val_A2_CHX_Input $max_val_A2_Puro $max_val_B1_CHX_mock $max_val_B1_CHX_Input $max_val_B1_Puro)
        # max=$(printf "%s\n" "${numbers[@]}" | sort -nr | head -n1)


        # max=20
        # echo $max
        # sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini
        # sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini
        # sed -i "s/^max_value *= *.*/max_value = $max/"$bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini

        # sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini
        # sed -i "s/^max_value *= *.*/max_value = $max/" $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini
        # sed -i "s/^max_value *= *.*/max_value = $max/"$bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini


    sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini
    pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini\
    --region $MANE_trans:0-$trans_length\
    --height 40\
    --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candidates_CHX_Impb_over_Mock_ratio.pdf


    sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini
    pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini\
    --region $MANE_trans:0-$trans_length\
      --height 40\
    --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candidates_Puro_Impb_over_Input_ratio.pdf


    sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini
    pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini\
        --region $MANE_trans:0-$trans_length\
        --height 40\
        --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candidates_CHX_Impb_over_Input_ratio.pdf



    # set number of bin to length of the transcript
    sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini
    pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini\
    --region $MANE_trans:0-$trans_length\
    --height 40\
    --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candiates_CHX_Impa_over_Mock_ratio.pdf

    
    sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini
    pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini\
    --region $MANE_trans:0-$trans_length\
    --height 40\
    --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candiates_Puro_Impa_over_Input_ratio.pdf

    sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini
    pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini\
        --region $MANE_trans:0-$trans_length\
        --height 40\
        --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candiates_CHX_Impa_over_Input_ratio.pdf
done