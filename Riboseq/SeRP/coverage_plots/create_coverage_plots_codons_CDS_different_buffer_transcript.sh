#!/bin/bash

bam_dir=$1
outdir=$2
coverage_script_dir=$3
mane_gtf=$4


if [ ! -d "$bam_dir"/enrichment_plots_CDS ]; then
         mkdir $bam_dir/enrichment_plots_CDS
fi

if [ ! -d "$outdir" ]; then
         mkdir $outdir
fi




################################################################################
# Filter BAMS for CDS coordinates                                              #
################################################################################
if [[ ! -e "$outdir/MANE_transcripts_114.txt" ]]; then
        cp /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots_CDS/CDS_coordinates/MANE_transcripts_114.txt "$outdir/MANE_transcripts_114.txt"
fi 
# # Do not merge replicates, but average them later one
# # single files are all sorted already
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

if [ ! -d "$bam_dir"/enrichment_plots_CDS/whole_transcript_bigwig ]; then
        mkdir $bam_dir/enrichment_plots_CDS/whole_transcript_bigwig

        bash ${coverage_script_dir}/BamCompare_Imp_vs_background.sh \
        $bam_dir \
        "CHX" \
        "In" \
        ""



        bash ${coverage_script_dir}/BamCompare_Imp_vs_background.sh \
        $bam_dir \
        "CHX" \
        "M" \
        "E"



        bash ${coverage_script_dir}/BamCompare_Imp_vs_background.sh \
        $bam_dir \
        "Puro" \
        "In" \
        ""

# Mock over Input need to do separately as not included in the sytax
        export TMPDIR=/scratch/tmp/$USER
        for bam in $bam_dir/*.bam; do
        mock_filename=$(basename "$bam")
        if [[ "$mock_filename" =~ ^uf_muellermcnicoll_([0-9_]+)_RR_M_CHX_E([0-9]+)\.cut\.fastp\.bowtie1_concat_transcriptome_k1_R1_sorted_filtered_q10\.bam$ ]]; then
                date="${BASH_REMATCH[1]}"
                batch="${BASH_REMATCH[2]}"


                input_filename=$bam_dir/uf_muellermcnicoll_2025_05_??_RR_In_CHX_${batch}.cut.fastp.bowtie1_concat_transcriptome_k1_R1_sorted_filtered_q10.bam
                echo $input_filename
                echo $mock_filename
                
                bamCompare -b1 $bam_dir/$mock_filename\
                -b2 $input_filename\
                -o $bam_dir/enrichment_plots_CDS/whole_transcript_bigwig/M_E${batch}_over_In_CHX_whole_trans_b1_no_smooth.bw\
                --operation ratio \
                --binSize 1 \
                -of bigwig\
                -p 64

        fi
        done
fi




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
        "${bam_dir}/enrichment_plots_CDS/impA_0_5_M_and_input_DEG_plots"
        "${bam_dir}/enrichment_plots_CDS/impB_0_5_M_and_input_DEG_plots"
        "${bam_dir}/enrichment_plots_CDS/RIP_LFC2_IP"
        "${bam_dir}/enrichment_plots_CDS/SeRP_KPNAs"
        "${bam_dir}/enrichment_plots_CDS/impA_B_0_5_M_and_input_DEG_plots"
)


if [ ! -d "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig/coordinates_per_transcript_csvs" ]; then
        mkdir "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig/coordinates_per_transcript_csvs"
fi

# Loop through all
for i in "${!transcripts_to_plot_txts[@]}"; do
     to_plot_txt=${transcripts_to_plot_txts[$i]}  
     outdir_plot=${out_dirs_for_plot_files[$i]}

     if [ ! -d $outdir_plot ]; then
        mkdir $outdir_plot

        ############## Buffer 1's ####################################################
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "A2" \
        "In" \
        --puro "" \
        --color "#1eb0e6" \
        --window_size 63

        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "B1" \
        "In" \
        --puro "" \
        --color "#29449c" \
        --window_size 63

        # A2 over Mock
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "A2" \
        "M" \
        --puro "" \
        --color "#1eb0e6" \
        --window_size 63

        # B1 over Mock
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "B1" \
        "M" \
        --puro "" \
        --color "#29449c" \
        --window_size 63



        # Puro A2 over In
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "A2" \
        "In" \
        --puro "Puro" \
        --color "#e53e3e" \
        --window_size 63

        # Puro B1 over In
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "B1" \
        "In" \
        --puro "Puro" \
        --color "#9b2020" \
        --window_size 63


        # Mock over Input
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_one_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        ${to_plot_txt} \
        ${outdir_plot} \
        "M" \
        "In" \
        --puro "" \
        --color "#6d6d6d" \
        --window_size 63
     fi

     if [ ! -d ${outdir_plot}_buffer_UTR ]; then
        mkdir ${outdir_plot}_buffer_UTR

        ############# buffer UTR #####################################################
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_UTR_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_transcript_length.txt" \
        ${to_plot_txt} \
        ${outdir_plot}_buffer_UTR \
        "A2" \
        "In" \
        --puro "" \
        --color "#1eb0e6" \
        --window_size 63

        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_UTR_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_transcript_length.txt" \
        ${to_plot_txt} \
        ${outdir_plot}_buffer_UTR \
        "B1" \
        "In" \
        --puro "" \
        --color "#29449c" \
        --window_size 63

        # A2 over Mock
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_UTR_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_transcript_length.txt" \
        ${to_plot_txt} \
        ${outdir_plot}_buffer_UTR \
        "A2" \
        "M" \
        --puro "" \
        --color "#1eb0e6" \
        --window_size 63

        # B1 over Mock
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_UTR_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_transcript_length.txt" \
        ${to_plot_txt} \
        ${outdir_plot}_buffer_UTR \
        "B1" \
        "M" \
        --puro "" \
        --color "#29449c" \
        --window_size 63



        # Puro A2 over In
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_UTR_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_transcript_length.txt" \
        ${to_plot_txt} \
        ${outdir_plot}_buffer_UTR \
        "A2" \
        "In" \
        --puro "Puro" \
        --color "#e53e3e" \
        --window_size 63

        # Puro B1 over In
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_UTR_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_transcript_length.txt" \
        ${to_plot_txt} \
        ${outdir_plot}_buffer_UTR \
        "B1" \
        "In" \
        --puro "Puro" \
        --color "#9b2020" \
        --window_size 63


        # Mock over Input
        python ${coverage_script_dir}/pyBigWig_for_plotting_with_errors_UTR_buffer.py \
        "${bam_dir}/enrichment_plots_CDS/whole_transcript_bigwig" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_CDS_coordinates.bed" \
        "${bam_dir}/enrichment_plots_CDS/CDS_coordinates/MANE_transcript_length.txt" \
        ${to_plot_txt} \
        ${outdir_plot}_buffer_UTR \
        "M" \
        "In" \
        --puro "" \
        --color "#6d6d6d" \
        --window_size 63
     fi


done


coordinate_dir="/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/bowtie1/filtered/q10/enrichment_plots_CDS/whole_transcript_bigwig/coordinates_per_transcript_csvs"
# this directory has every transcript several times, 1 buffering, UTR buffering and the different conditions A, B, M, Puro
# idea: just loop through all files: pd_df: file name: then it is clear what it belongs to, redundancy does not matter

if [ ! -d "${coordinate_dir}/onset_caclulation" ]; then
         mkdir ${coordinate_dir}/onset_caclulation
fi


# enrichment_length: 51:3, 17 AA or 51 bp
# change to 50 AA
python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 17 \
        --enrichment_threshold 2


python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 50 \
        --enrichment_threshold 2


python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 75 \
        --enrichment_threshold 2

python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 100 \
        --enrichment_threshold 2



# enrichment_length: 51:3, 17 AA or 51 bp
# change to 50 AA

python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 50 \
        --enrichment_threshold 1.5


python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 75 \
        --enrichment_threshold 1.5

python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 100 \
        --enrichment_threshold 1.5


# enrichment_length: 51:3, 17 AA or 51 bp
# change to 50 AA

python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 50 \
        --enrichment_threshold 1.75


python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 75 \
        --enrichment_threshold 1.75

python ${coverage_script_dir}/get_onset.py \
        ${coordinate_dir} \
        --enrichment_length 100 \
        --enrichment_threshold 1.75