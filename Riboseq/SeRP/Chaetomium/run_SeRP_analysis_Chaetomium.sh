#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Riboseq

script_dir="/home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1"

indir="/projects/serp/work/data/SeRP_April_2025/Chaetomium"
outdir="/projects/serp/work/Output/April_2025/Chaetomium"
outdir_FASTQC1=${outdir}"/preprocess/fastqc_unprocessed"
outdir_CUTADAPT=${outdir}"/preprocess/cutadapt"
fastpOut=${outdir}"/preprocess/fastp"
fastpFASTQC=${outdir}"/preprocess/fastp/fastqc"
# Bowtie2_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
# Bowtie2_base_name="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/index"
# Bowtie2_out_dir="/projects/serp/work/Output/April_2025/Chaetomium/transcriptome_mapping"
# Genome_Fasta=
# EnsemblFilteredRef=
bowtie_outdir="/projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome"

if [ ! -d $outdir/preprocess ]; then
        mkdir $outdir/preprocess
fi

if [ ! -d $outdir_FASTQC1 ]; then
        mkdir $outdir_FASTQC1
fi

if [ ! -d $outdir_CUTADAPT ]; then
        mkdir $outdir_CUTADAPT
fi


if [ ! -d $fastpOut ]; then
        mkdir $fastpOut
fi


if [ ! -d "out_reports_of_runs" ]; then
        mkdir out_reports_of_runs
fi

################################################################################
# QC and PREPROCESSING                                                         #
################################################################################

# bash ../preprocessing_cutadapt_steps_importins.sh \
#  $indir \
#   $outdir_FASTQC1 \
#   $outdir_CUTADAPT \
#   $fastpOut \
#   $fastpFASTQC\
# > "out_reports_of_runs/preprocessing_cutadapt_chaetomium.out" 2>&1

# python ${script_dir}/preprocessing/cutadapt_output_parsing.py \
#  "out_reports_of_runs/preprocessing_cutadapt_chaetomium.out" \
#   "/home/ckalk/scripts/SplitORFs/Riboseq/SeRP/out_reports_of_runs/cutadapt_summary.csv"




################################################################################
# TRANSCRIPTOMIC ALIGNMENT                                                     #
################################################################################
# filter transcriptome for the longest transcript per gene
# first sort how cgat requires it

conda activate pygtftk
python filter_gtf_for_longest_transcripts.py \
 /projects/serp/work/references/Supplementary_File_2.gtf \
 /projects/serp/work/references/Supplementary_File_2_longest_transcript.gtf


conda activate Riboseq
# generate transcriptome
gffread /projects/serp/work/references/Supplementary_File_2_longest_transcript.gtf \
 -g /projects/serp/work/references/Chaetomium_thermophilum_var_thermophilum_dsm_1495.CTHT_3.0.dna.toplevel.fa \
 -w /projects/serp/work/references/Chaetomium_thermophilum_longest_transcript.fasta



# align to transcriptome
bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/bowtie1_align_21_10_25.sh \
 ${bowtie_outdir}/index \
 /projects/serp/work/references/Chaetomium_thermophilum_longest_transcript.fasta \
 ${fastpOut} \
 ${bowtie_outdir} \
 Chaetomium_transcriptome \
 /home/ckalk/scripts/SplitORFs/Riboseq/SeRP/Chaetomium/out_reports_of_runs/Chaetomium_align_bowtie.out \
 > /home/ckalk/scripts/SplitORFs/Riboseq/SeRP/Chaetomium/out_reports_of_runs/Chaetomium_align_bowtie.out 2>&1


python get_trna_rrna_tids.py

python summarize_bowtie2_alns_by_source_chaetomium_noncoding.py \
    /projects/serp/work/references/Chaetomium_thermophilum_noncoding.fasta,/projects/serp/work/references/Chaetomium_thermophilum_protein_coding.fasta \
    ${bowtie_outdir}

python summarize_bowtie2_alns_by_source_chaetomium_noncoding.py \
    /projects/serp/work/references/Chaetomium_thermophilum_noncoding.fasta,/projects/serp/work/references/Chaetomium_thermophilum_protein_coding.fasta \
    ${bowtie_outdir}/filtered

python summarize_bowtie2_alns_by_source_chaetomium_noncoding.py \
    /projects/serp/work/references/Chaetomium_thermophilum_noncoding.fasta,/projects/serp/work/references/Chaetomium_thermophilum_protein_coding.fasta \
    ${bowtie_outdir}/filtered/q10


if [ ! -d ${bowtie_outdir}/filtered/q10/DEGs ]; then
        mkdir ${bowtie_outdir}/filtered/q10/DEGs
fi

samtools faidx "/projects/serp/work/references/Chaetomium_thermophilum_protein_coding.fasta"
Rscript  PCA_conditions_DeSeq2_SeRP_Chaetomium.R \
"${bowtie_outdir}"/filtered/q10


 grep -Fxf ${bowtie_outdir}/filtered/q10/DEGs/E_S_over_E_WT_0_5_0.05.txt \
 ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_0_5_0.05.txt \
> ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_0_5_and_E_S_over_E_WT_0_5.txt


 grep -Fxf ${bowtie_outdir}/filtered/q10/DEGs/E_S_over_E_WT_1_0_0.01.txt \
 ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_1_0_0.01.txt \
> ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0.txt

# -v invert match only select lines that do no match
 grep -Fxvf  ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_WT_0_5_0.05.txt \
 ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_0_5_and_E_S_over_E_WT_0_5.txt \
> ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_0_5_and_E_S_over_E_WT_0_5_not_IP_WT_vs_IN.txt


 grep -Fxvf ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_WT_1_0_0.01.txt \
 ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0.txt \
> ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0_not_IP_WT_vs_IN.txt


 grep -Fxvf ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_WT_0_5_0.05.txt \
 ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0.txt \
> ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0_not_IP_WT_vs_IN_0.5_andpadj_0.05.txt

python compare_old_new_results.py \
        /projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/Analysis_31_10_25_Remus/final_results_SND3_WT_LFC1_padj0.01.csv  \
        /projects/serp/work/Output/April_2025/Chaetomium/align_transcriptome/filtered/q10/DEGs/Analysis_31_10_25_Remus/final_results_SND3_WT_LFC1_padj0.01_all.csv \
        ${bowtie_outdir}/filtered/q10/DEGs/E_over_In_S_1_0_and_E_S_over_E_WT_1_0_not_IP_WT_vs_IN.txt

# conda activate Riboseq
# if [ ! -d "${bowtie_outdir}"/filtered/q10/Ribowaltz ]; then
#         "${bowtie_outdir}"/filtered/q10/Ribowaltz
# fi
# Rscript RiboWaltz_SeRP_Chaetomium_single_samples.R \
# "${bowtie_outdir}"/filtered/q10



################################################################################
# SeRP coverage plots                                                          #
################################################################################
# bash create_coverage_plots_importins_over_mock.sh \
#  "${Bowtie2_out_dir}"/filtered/q10

# bash ${coverage_script_dir}/create_coverage_plots_codons_whole_transcript.sh \
#     "${Bowtie2_out_dir}"/filtered/q10 \
#     ""${Bowtie2_out_dir}"/filtered/q10/enrichment_plots_CDS/CDS_coordinates" \
#     ${coverage_script_dir} \
#     $mane_gtf

################################################################################
# Check al counts with bam multicov                                            #
################################################################################

# faidx --transform bed ${Bowtie2_ref_fasta} > $(dirname ${Bowtie2_ref_fasta})/Ignolia_transcriptome_and_contamination.bed

# reference_bed=$(dirname ${Bowtie2_ref_fasta})/Ignolia_transcriptome_and_contamination.bed

# for bam in "${Bowtie2_out_dir}"/filtered/q10/*.bam; do
    # bedtools multicov -bams $bam -bed $reference_bed > "${Bowtie2_out_dir}"/filtered/q10/$(basename $bam .bam)_mulitcov.bed
#     diff <(cut -f 1,4 "${Bowtie2_out_dir}"/filtered/q10/$(basename $bam .bam)_mulitcov.bed)\
#      <(cut -f 1,3  "${Bowtie2_out_dir}"/filtered/q10/$(basename $bam .bam)_idxstats.out)
# done














