#!/bin/bash

bam_dir=$1
outdir=$2
coverage_script_dir=$3
mane_gtf=$4


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
# ${mane_gtf} \
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

# Mock over Input need to do separately as not included in the sytax
# export TMPDIR=/scratch/tmp/$USER
# for bam in $bam_dir/*.bam; do
#   mock_filename=$(basename "$bam")
#   if [[ "$mock_filename" =~ ^uf_muellermcnicoll_([0-9_]+)_RR_M_CHX_E([0-9]+)\.cut\.fastp\.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10\.bam$ ]]; then
#         date="${BASH_REMATCH[1]}"
#         batch="${BASH_REMATCH[2]}"


#         input_filename=$bam_dir/uf_muellermcnicoll_2025_05_??_RR_In_CHX_${batch}.cut.fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10.bam
#         echo $input_filename
#         echo $mock_filename
        
#         bamCompare -b1 $bam_dir/$mock_filename\
#          -b2 $input_filename\
#          -o $bam_dir/enrichment_plots_CDS/whole_transcript_bigwig/M_E${batch}_over_In_CHX_whole_trans_b1_no_smooth.bw\
#          --operation ratio \
#          --binSize 1 \
#          -of bigwig\
#          -p 64

#   fi
# done




# for pyBigWig plotting script
CDS_coordinates="$outdir/MANE_CDS_coordinates.bed"


# python map_gids_to_MANE_tids.py \
#  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
#  $bam_dir/enrichment_plots/2025-06-04_RIP-Seq_hits.csv  \
#  $bam_dir/enrichment_plots/RIP_hits_MANE_tIDs.txt \
#  $bam_dir/enrichment_plots/RIP_hits_gids_to_MANE_tIDs.csv \
#  multirow



################################################################################
# Imp a enrichment plots                                                       #
################################################################################

transcripts_to_plot_txts=(
        "${bam_dir}/DEGs/DEGs_A2_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs.txt"
        "${bam_dir}/DEGs/DEGs_B1_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs.txt"
        "${bam_dir}/DEGs/RIP_LFC2_IP_against_Input_MANE_tIDs.txt"
        "${bam_dir}/DEGs/SeRP_KPNAs_MANE_tIDs_modified.txt"
        "${bam_dir}/DEGs/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt"
        )

out_dirs_for_plot_files=(
        "${bam_dir}/enrichment_plots_CDS/impA_0_5_M_and_input_DEG_plots_whole_transcript"
        "${bam_dir}/enrichment_plots_CDS/impB_0_5_M_and_input_DEG_plots_whole_transcript"
        "${bam_dir}/enrichment_plots_CDS/RIP_LFC2_IP"
        "${bam_dir}/enrichment_plots_CDS/SeRP_KPNAs"
        "${bam_dir}/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript"
)

# Loop through all
for i in "${!transcripts_to_plot_txts[@]}"; do
     to_plot_txt=${transcripts_to_plot_txts[$i]}  
     outdir_plot=${out_dirs_for_plot_files[$i]}

     if [ ! -d $outdir_plot ]; then
        mkdir $outdir_plot
     fi

     if [ ! -d $outdir_plot/plot_input_and_mock ]; then
        mkdir $outdir_plot/plot_input_and_mock
     fi

        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "A2" \
        "In" \
        --puro "" \
        --color "#1eb0e6"

        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "B1" \
        "In" \
        --puro "" \
        --color "#29449c"

        # python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript_min_max.py \
        # "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        # "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        # ${to_plot_txt} \
        # ${outdir_plot} \
        # "A2" \
        # "In" \
        # --puro "" \
        # --color "#1eb0e6"

        # python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript_min_max.py \
        # "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        # "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        # ${to_plot_txt} \
        # ${outdir_plot} \
        # "B1" \
        # "In" \
        # --puro "" \
        # --color "#29449c"



        # A2 over Mock
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "A2" \
        "M" \
        --puro "" \
        --color "#1eb0e6"

        # B1 over Mock
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "B1" \
        "M" \
        --puro "" \
        --color "#29449c"

        # Puro A2 over In
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "A2" \
        "In" \
        --puro "Puro" \
        --color "#1eb0e6"

        # Puro B1 over In
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "B1" \
        "In" \
        --puro "Puro" \
        --color "#29449c"


        # Mock over Input
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_whole_transcript.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "M" \
        "In" \
        --puro "" \
        --color "#6d6d6d"


        readarray -t trans_to_plot < ${to_plot_txt}

        for MANE_trans in "${trans_to_plot[@]}"
        do
                if grep -q -F "${MANE_trans}" "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed"; then
                        python ${coverage_script_dir}/plot_two_comparisons_in_one_plot.py \
                        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig/coordinates_per_transcript_csvs" \
                        $MANE_trans \
                        $outdir_plot/plot_input_and_mock \
                        'A2' \
                        --background1 'In' \
                        --background2 'M' \
                        --puro '' \
                        --color1 '#1eb0e6' \
                        --color2 'orange'
                fi
        done
done


# readarray -t importin_enriched_A2_and_B1 < ${outdir}/DEGs/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt

# for MANE_trans in "${importin_enriched_A2_and_B1[@]}"
# do
#         python ${coverage_script_dir}/plot_two_comparisons_in_one_plot.py \
#         '${outdir}/enrichment_plots_CDS/whole_transcript_bigwig/coordinates_per_transcript_csvs' \
#         $MANE_trans \
#         '${outdir}/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots_whole_transcript/plot_input_and_mock' \
#         'A2' \
#         --background1 'In' \
#         --background2 'M' \
#         --puro '' \
#         --color1 '#1eb0e6' \
#         --color2 'orange'

# done