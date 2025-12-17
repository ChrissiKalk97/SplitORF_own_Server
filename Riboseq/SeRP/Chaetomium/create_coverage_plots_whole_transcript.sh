#!/bin/bash

bam_dir=$1
outdir=$2
coverage_script_dir=$3


################################################################################
# create bigwig whole transcript                                               #
################################################################################
eval "$(conda shell.bash hook)"
conda activate Riboseq
export TMPDIR=/scratch/tmp/$USER
if [ ! -d "$outdir" ]; then
         mkdir $outdir
    # E S over In S
    bash ${coverage_script_dir}/BamCompare_Chaetomium.sh \
        $bam_dir \
        In \
        E \
        S \
        S \
        $outdir



    # E S over E WT
    bash ${coverage_script_dir}/BamCompare_Chaetomium.sh \
        $bam_dir \
        E \
        E \
        W \
        S \
        $outdir



    # E WT over In WT
    bash ${coverage_script_dir}/BamCompare_Chaetomium.sh \
        $bam_dir \
        In \
        E \
        W \
        W \
        $outdir
fi


################################################################################
# S E over S In enrichment plots                                               #
################################################################################

transcripts_to_plot_txts=(
        "${bam_dir}/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0_not_IP_WT_vs_IN_0.5_andpadj_0.05.txt"
        )

out_dirs_for_plot_files=(
        "${outdir}/E_over_In_S_1_0_and_E_S_over_E_WT_1_0_not_IP_WT_vs_IN_0.5_DEG_plots"
)

if [ ! -d $outdir/coordinates_per_transcript_csvs ]; then
    mkdir $outdir/coordinates_per_transcript_csvs
fi

# Loop through all
for i in "${!transcripts_to_plot_txts[@]}"; do
     to_plot_txt=${transcripts_to_plot_txts[$i]}  
     outdir_plot=${out_dirs_for_plot_files[$i]}

        if [ ! -d ${outdir_plot} ]; then
            mkdir $outdir_plot
        fi 

        # path_to_bw_files = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/enrichment_plots_whole_trans'
        # transcript_fai = '/projects/serp/work/references/Chaetomium_thermophilum_longest_transcript.fasta.fai'
        # transcripts_to_plot_txt = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0_not_IP_WT_vs_IN_0.5_andpadj_0.05.txt'
        # out_path = '/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/enrichment_plots_whole_trans/S_E_over_S_In_DEG_plots'
        # numerator = 'S_E'
        # background = 'S_In'
        # color = '#1eb0e6'
        # window_size = 63
        ############## Buffer 1's ####################################################
        # python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer_Chaetomium.py \
        # "${bam_dir}/enrichment_plots_whole_trans" \
        # "/projects/serp/work/references/Chaetomium_thermophilum_longest_transcript.fasta.fai" \
        # ${to_plot_txt} \
        # ${outdir_plot} \
        # --numerator "S_E" \
        # --background "S_In" \
        # --color "#1eb0e6" \
        # --window_size 63

        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer_Chaetomium.py \
        "${bam_dir}/enrichment_plots_whole_trans" \
        "/projects/serp/work/references/Chaetomium_thermophilum_longest_transcript.fasta.fai" \
        ${to_plot_txt} \
        ${outdir_plot} \
        --numerator "S_E" \
        --background "W_E" \
        --color "#1eb0e6" \
        --window_size 63

        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer_Chaetomium.py \
        "${bam_dir}/enrichment_plots_whole_trans" \
        "/projects/serp/work/references/Chaetomium_thermophilum_longest_transcript.fasta.fai" \
        ${to_plot_txt} \
        ${outdir_plot} \
        --numerator "W_E" \
        --background "W_In" \
        --color "#1eb0e6" \
        --window_size 63

done




# coordinate_dir="/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/bowtie1/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig/coordinates_per_transcript_csvs"
# # this directory has every transcript several times, 1 buffering, UTR buffering and the different conditions A, B, M, Puro
# # idea: just loop through all files: pd_df: file name: then it is clear what it belongs to, redundancy does not matter

# if [ ! -d "${coordinate_dir}/onset_caclulation" ]; then
#          mkdir ${coordinate_dir}/onset_caclulation
# fi


# # enrichment_length: 51:3, 17 AA or 51 bp
# # change to 50 AA
# python ${coverage_script_dir}/get_onset.py \
#         ${coordinate_dir} \
#         --enrichment_length 17 \
#         --enrichment_threshold 2


# python ${coverage_script_dir}/get_onset.py \
#         ${coordinate_dir} \
#         --enrichment_length 50 \
#         --enrichment_threshold 2


# python ${coverage_script_dir}/get_onset.py \
#         ${coordinate_dir} \
#         --enrichment_length 75 \
#         --enrichment_threshold 2

# python ${coverage_script_dir}/get_onset.py \
#         ${coordinate_dir} \
#         --enrichment_length 100 \
#         --enrichment_threshold 2