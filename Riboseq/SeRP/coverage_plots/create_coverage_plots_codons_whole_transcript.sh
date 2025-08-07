#!/bin/bash

bam_dir=$1
outdir=$2
coverage_script_dir=$3

# need to supply these via the main SeRP call
# '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig'
#  '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed'
# '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt'
# '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript' 
# "/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf"

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
# python ${coverage_script_dir}/filter_Ensembl_structures_for_CDS_coords.py \
# "$outdir/MANE_transcripts_114.txt" \
# "/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf" \
# "$outdir/MANE_CDS_coordinates"



################################################################################
# create bigwig whole transcript                                               #
################################################################################
eval "$(conda shell.bash hook)"
conda activate Riboseq
export TMPDIR=/scratch/tmp/$USER

# bash ${coverage_script_dir}/BamCompare_Imp_vs_background.sh \
#  $bam_dir \
#  "CHX" \
#  "In" \
#  ""



#  bash ${coverage_script_dir}/BamCompare_Imp_vs_background.sh \
#  $bam_dir \
#  "CHX" \
#  "M" \
#  "E"



#  bash ${coverage_script_dir}/BamCompare_Imp_vs_background.sh \
#  $bam_dir \
#  "Puro" \
#  "In" \
#  ""




# for pyBigWig plotting script
CDS_coordinates="$outdir/MANE_CDS_coordinates.bed"


# python map_gids_to_MANE_tids.py \
#  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
#  $bam_dir/enrichment_plots/2025-06-04_RIP-Seq_hits.csv  \
#  $bam_dir/enrichment_plots/RIP_hits_MANE_tIDs.txt \
#  $bam_dir/enrichment_plots/RIP_hits_gids_to_MANE_tIDs.csv \
#  multirow





if [ ! -d $bam_dir/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots ]; then
        mkdir $bam_dir/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots
fi

if [ ! -d $bam_dir/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript ]; then
        mkdir $bam_dir/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript
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


python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript.py \
 '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig' \
 '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed' \
'/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt' \
'/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript' \
'A2' \
 'In' \
 --puro '' \
 --color '#1eb0e6'

python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript.py \
 '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig' \
 '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed' \
'/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt' \
'/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript' \
'B1' \
 'In' \
 --puro '' \
 --color '#29449c'

 python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript_min_max.py \
 '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig' \
 '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed' \
'/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt' \
'/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript' \
'A2' \
 'In' \
 --puro '' \
 --color '#1eb0e6'

python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript_min_max.py \
 '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig' \
 '/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed' \
'/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt' \
'/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript' \
'B1' \
 'In' \
 --puro '' \
 --color '#29449c'



################################################################################
# Imp b enrichment plots                                                       #
################################################################################

# readarray -t importin_enriched_B1 < /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_B1_CHX_enriched_over_Input_and_Mock_MANE_tIDs.txt


################################################################################
# plot hot candidates on the same scale                                        #
################################################################################
# hot_candidates=(ENST00000397885.3 ENST00000330560.8 ENST00000428849.7 ENST00000611405.5)



################################################################################
# plot histones and GAPDH                                                      #
################################################################################
# readarray -t histone_array < $bam_dir/enrichment_plots/RIP_hits_MANE_tIDs_GAPDH.txt

