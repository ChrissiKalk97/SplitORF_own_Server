#!/bin/bash

bam_dir=$1
outdir=$2


if [ ! -d $bam_dir/enrichment_plots_CDS ]; then
         mkdir $bam_dir/enrichment_plots_CDS
fi

if [ ! -d $outdir ]; then
         mkdir $outdir
fi


################################################################################
# Filter BAMS for CDS coordinates                                              #
################################################################################
# Do not merge replicates, but average them later one
# single files are all sorted already
# eval "$(conda shell.bash hook)"
# conda activate pygtftk
# python /home/ckalk/scripts/SplitORFs/Riboseq/SeRP/coverage_plots/filter_Ensembl_structures_for_CDS_coords.py \
# "$outdir/MANE_transcripts_114.txt" \
# "/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf" \
# "$outdir/MANE_CDS_coordinates"

CDS_coordinates="$outdir/MANE_CDS_coordinates.bed"

if [ ! -d "${outdir}"/CDS_bams/ ]; then
        mkdir $outdir/CDS_bams/
fi

eval "$(conda shell.bash hook)"
conda activate Riboseq
# for bam in $bam_dir/*.bam; do
#         sample=$(basename $bam .cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10.bam)
#         bedtools intersect -a $bam -b $CDS_coordinates > $outdir/CDS_bams/${sample}_CDS.bam
#         samtools sort $outdir/CDS_bams/${sample}_CDS.bam > $outdir/CDS_bams/${sample}_CDS_sorted.bam
#         samtools index $outdir/CDS_bams/${sample}_CDS_sorted.bam
# done





################################################################################
# create CHX over Input bigwig                                                 #
################################################################################
export TMPDIR=/scratch/tmp/$USER
for bam in $outdir/CDS_bams/*.bam; do
  importin_filename=$(basename "$bam")
  if [[ "$importin_filename" =~ ^uf_muellermcnicoll_([0-9_]+)_RR_([AB][12])_CHX_E([0-9]+)_CDS_sorted\.bam$ ]]; then
        date="${BASH_REMATCH[1]}"
        importin="${BASH_REMATCH[2]}"
        batch="${BASH_REMATCH[3]}"


        input_filename=$outdir/CDS_bams/uf_muellermcnicoll_2025_05_??_RR_In_CHX_"$batch"_CDS_sorted.bam
        echo $input_filename
        echo $importin_filename
        
        bamCompare -b1 $outdir/CDS_bams/$importin_filename\
         -b2 $input_filename\
         -o $bam_dir/enrichment_plots_CDS/${importin}_E${batch}_over_Input_CDS_b1_no_smooth.bw\
         --operation ratio \
         --binSize 1 \
         #--smoothLength 21 \
         -of bigwig\
         -p 64

  fi
done
wait

################################################################################
# create CHX over Mock bigwig                                                 #
################################################################################
# for bam in $outdir/CDS_bams/*.bam; do
#   importin_filename=$(basename "$bam")
#   if [[ "$importin_filename" =~ ^uf_muellermcnicoll_([0-9_]+)_RR_([AB][12])_CHX_E([0-9]+)_CDS_sorted\.bam$ ]]; then
#         date="${BASH_REMATCH[1]}"
#         importin="${BASH_REMATCH[2]}"
#         batch="${BASH_REMATCH[3]}"


#         mock_filename=$outdir/CDS_bams/uf_muellermcnicoll_2025_05_??_RR_M_CHX_E"$batch"_CDS_sorted.bam
#         echo $input_filename
#         echo $mock_filename
        
#         bamCompare -b1 $outdir/CDS_bams/$importin_filename\
#          -b2 $mock_filename\
#          -o $bam_dir/enrichment_plots_CDS/${importin}_E${batch}_over_Mock_CDS_b3_sl_21.bw\
#          --operation ratio \
#          --binSize 3 \
#          --smoothLength 21 \
#          -of bigwig\
#          -p 32
#   fi
# done
# wait 


################################################################################
# create Puro over Input bigwig                                                 #
################################################################################
# for bam in $outdir/CDS_bams/*.bam; do
#   {
#   importin_filename=$(basename "$bam")
#   if [[ "$importin_filename" =~ ^uf_muellermcnicoll_([0-9_]+)_RR_([AB][12])_CHX_E([0-9]+)_CDS_sorted\.bam$ ]]; then
#         date="${BASH_REMATCH[1]}"
#         importin="${BASH_REMATCH[2]}"
#         batch="${BASH_REMATCH[3]}"


#         mock_filename=$outdir/CDS_bams/uf_muellermcnicoll_2025_05_??_RR_M_CHX_"$batch"_CDS_sorted.bam
#         echo $input_filename
#         echo $importin_filename
        
#         bamCompare -b1 $outdir/CDS_bams/$importin_filename\
#          -b2 $mock_filename\
#          -o $bam_dir/enrichment_plots_CDS/${importin}_E${batch}_over_Mock_CDS_b3_sl_21.bw\
#          --operation ratio \
#          --binSize 3 \
#          --smoothLength 21 \
#          -of bigwig\
#          -p 4
#   }&
#   fi
# done
# wait 




# python map_gids_to_MANE_tids.py \
#  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
#  $bam_dir/enrichment_plots/2025-06-04_RIP-Seq_hits.csv  \
#  $bam_dir/enrichment_plots/RIP_hits_MANE_tIDs.txt \
#  $bam_dir/enrichment_plots/RIP_hits_gids_to_MANE_tIDs.csv \
#  multirow





if [ ! -d $bam_dir/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots ]; then
        mkdir $bam_dir/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots
fi

# if [ ! -d $bam_dir/enrichment_plots/impB_DEG_plots ]; then
#         mkdir $bam_dir/enrichment_plots/impB_DEG_plots
# fi

# if [ ! -d $bam_dir/enrichment_plots/hot_candiates ]; then
#         mkdir $bam_dir/enrichment_plots/hot_candiates
# fi

# if [ ! -d $bam_dir/enrichment_plots/histones ]; then
#         mkdir $bam_dir/enrichment_plots/histones
# fi


################################################################################
# Imp a enrichment plots                                                       #
################################################################################

readarray -t importin_enriched_A2 < /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt


if [ ! -d $bam_dir/enrichment_plots_CDS/bwtool_output ]; then
        mkdir $bam_dir/enrichment_plots_CDS/bwtool_output
fi

# get the raw tab file
# multiBigwigSummary bins -b $bam_dir/enrichment_plots_CDS/A2_E1_over_Input_CDS_b3_sl_21.bw \
# $bam_dir/enrichment_plots_CDS/A2_E2_over_Input_CDS_b3_sl_21.bw \
# $bam_dir/enrichment_plots_CDS/A2_E3_over_Input_CDS_b3_sl_21.bw \
# --binSize 3 \
# -p 16 \
# -o $bam_dir/enrichment_plots_CDS/A2_average_over_Input_CDS_b3_sl_21.npz \
# --outRawCounts $bam_dir/enrichment_plots_CDS/A2_average_over_Input_CDS_b3_sl_21_raw.tab

# for MANE_trans in "${importin_enriched_A2[@]}"
# do
#         # python script: read tab file with pandas
#         # filter for MANE_transcript
#         # calcualte mean and std in pandas
#         # plot the 3 lines (mean, mean+std, mean-std)
#         # save figure

#         # try put bwtool
        

#       CDS_start=$(grep $MANE_trans $outdir/MANE_CDS_coordinates.bed | cut -f2)
#       CDS_end=$(grep $MANE_trans $outdir/MANE_CDS_coordinates.bed | cut -f3)
#       CDS_length=$((CDS_end - CDS_start))
#         echo $CDS_length

#         grep $MANE_trans $outdir/MANE_CDS_coordinates.bed > $outdir/${MANE_trans}_CDS.bed
#       /home/ckalk/bwtool_installation/bwtool/bwtool aggregate 0:$CDS_length $outdir/${MANE_trans}_CDS.bed \
#       $bam_dir/enrichment_plots_CDS/A2_E1_over_Input_CDS_b3_sl_21.bw,$bam_dir/enrichment_plots_CDS/A2_E2_over_Input_CDS_b3_sl_21.bw,$bam_dir/enrichment_plots_CDS/A2_E3_over_Input_CDS_b3_sl_21.bw \
#       $bam_dir/enrichment_plots_CDS/bwtool_output/${MANE_trans}_output.txt -starts -firstbase -expanded
#     # set number of bin to length of the transcript
# #     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini
# #     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini\
# #     --region $MANE_trans:0-$trans_length\
# #     --height 40\
# #     --outFileName $bam_dir/enrichment_plots/impA_DEG_plots/${MANE_trans}_CHX_Impa_over_Mock_ratio.pdf



        

    
# #     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini
# #     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini\
# #     --region $MANE_trans:0-$trans_length\
# #     --height 40\
# #     --outFileName $bam_dir/enrichment_plots/impA_DEG_plots/${MANE_trans}_Puro_Impa_over_Input_ratio.pdf

# #     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini
# #     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini\
# #         --region $MANE_trans:0-$trans_length\
# #         --height 40\
# #         --outFileName $bam_dir/enrichment_plots/impA_DEG_plots/${MANE_trans}_CHX_Impa_over_Input_ratio.pdf

# done


################################################################################
# Imp b enrichment plots                                                       #
################################################################################

# readarray -t importin_enriched_B1 < /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_B1_CHX_enriched_over_Input_and_Mock_MANE_tIDs.txt


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
# hot_candidates=(ENST00000397885.3 ENST00000330560.8 ENST00000428849.7 ENST00000611405.5)

# for MANE_trans in "${hot_candidates[@]}"
# do
#         trans_length=$(grep $MANE_trans  /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai | cut -f2)


#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candidates_CHX_Impb_over_Mock_ratio.pdf


#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#       --height 40\
#     --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candidates_Puro_Impb_over_Input_ratio.pdf


#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini\
#         --region $MANE_trans:0-$trans_length\
#         --height 40\
#         --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candidates_CHX_Impb_over_Input_ratio.pdf



#     # set number of bin to length of the transcript
#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candiates_CHX_Impa_over_Mock_ratio.pdf

    
#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candiates_Puro_Impa_over_Input_ratio.pdf

#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini\
#         --region $MANE_trans:0-$trans_length\
#         --height 40\
#         --outFileName $bam_dir/enrichment_plots/hot_candiates/${MANE_trans}_hot_candiates_CHX_Impa_over_Input_ratio.pdf
# done




################################################################################
# plot histones and GAPDH                                                      #
################################################################################
# readarray -t histone_array < $bam_dir/enrichment_plots/RIP_hits_MANE_tIDs_GAPDH.txt

# for MANE_trans in "${histone_array[@]}"
# do
#         trans_length=$(grep $MANE_trans  /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta.fai | cut -f2)


#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_CHX_over_mock_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/histones/${MANE_trans}_histones_CHX_Impb_over_Mock_ratio.pdf


#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_puro_over_in_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#       --height 40\
#     --outFileName $bam_dir/enrichment_plots/histones/${MANE_trans}_histones_Puro_Impb_over_Input_ratio.pdf


#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impb_CHX_over_in_ratio.ini\
#         --region $MANE_trans:0-$trans_length\
#         --height 40\
#         --outFileName $bam_dir/enrichment_plots/histones/${MANE_trans}_histones_CHX_Impb_over_Input_ratio.pdf



#     # set number of bin to length of the transcript
#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_mock_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/histones/${MANE_trans}_histones_CHX_Impa_over_Mock_ratio.pdf

    
#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_puro_over_in_ratio.ini\
#     --region $MANE_trans:0-$trans_length\
#     --height 40\
#     --outFileName $bam_dir/enrichment_plots/histones/${MANE_trans}_histones_Puro_Impa_over_Input_ratio.pdf

#     sed -i "s/^number_of_bins *= *.*/number_of_bins = $trans_length/" $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini
#     pyGenomeTracks --tracks $bam_dir/enrichment_plots/tracks_impa_CHX_over_in_ratio.ini\
#         --region $MANE_trans:0-$trans_length\
#         --height 40\
#         --outFileName $bam_dir/enrichment_plots/histones/${MANE_trans}_histones_CHX_Impa_over_Input_ratio.pdf
# done