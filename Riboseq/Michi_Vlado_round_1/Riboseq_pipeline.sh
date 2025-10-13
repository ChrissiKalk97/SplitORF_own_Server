#!/bin/bash


# =============================================================================
# Script Name: Riboseq_pipeline.sh
# Description: This pipeline contains all the steps for the analysis of  
#               different types of Riboseq data, all of the steps can 
#               chosen one by one or all at once:
#              - Step 1: Data preprocessing
#              - Step 2: Transcriptomic alignment (Ingolia reference)
#              - Step 3: Transcriptomic deduplication
#              - Step 4: Genomic alignment, deduplication
#              - Step 5: Split-ORF analysis
#               options: b: start from BAM files -> FASTQ
#                        d: deuplicate
#                        c: cutadapt
#                        s: soft-clip
#                        t: transcriptomic alignment
#                        g: genomic alignment
#               The pipeline needs to be run once per dataset
# Usage:       bash Riboseq_pipeline.sh args...
# Author:      Christina Kalk
# Date:        2025-10-10
# =============================================================================



eval "$(conda shell.bash hook)"
conda activate Riboseq

# Help message:
usage="
Usage: ./Riboseq_pipeline.sh [-options] [arguments]

where:
-h			show this help
"

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages


#available options for the programm
while getopts 'i:o:g:hbcdqr:stp' option; do
  case "$option" in
    i)
        indir="$OPTARG"
        ;;
    o)
        outdir="$OPTARG"
        ;;
    h) 
        echo "$usage"
        exit 1
        ;;
    b)
        start_bam=true
      ;;
    c) 
        cutadapt=true
      ;;
    d)
        dedup=true
      ;;
    g)
        genomic=true
        gtf="$OPTARG"
      ;;
    q)
        riboseqc=true
      ;;
    r)
        report_dir="$OPTARG"
      ;;
    s)
        soft_clip=true
      ;;
    t)
        transcriptomic=true
      ;;
    p)
        paired_reads=true
      ;;

   \?) 
        printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
done

# --- Check required options ---
if [[ -z "$indir" || -z "$outdir" ]]; then
  echo "Error: Both -i and -o options are required with arguments" >&2
  exit 1
fi

echo "All required options provided: -i=$indir -o=$outdir"


################################################################################
# PATH DEFINTIONS                                                              #
################################################################################

genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
bowtie2_ref_fasta="/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta"
bowtie2_base_name="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/index"

# bowtie2_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"


# output_star_transcriptomic="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"
# umi_dedup_outdir_transcriptomic="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/deduplicated"


################################################################################
# QC and PREPROCESSING                                                         #
################################################################################

if [[ $start_bam == true ]]; then
    bash /home/ckalk/scripts/SplitORFs/UPF1_deletion/data_download/get_fastq.sh $indir
    indir="${indir}/fastq"
fi

if [[ $cutadapt == true && $dedup == true && $paired_reads == true ]]; then
    outdir_preprocess=${outdir}/preprocess
    outdir_fastqc1=${outdir_preprocess}/fastqc
    outdir_cutadapt=${outdir_preprocess}/cutadapt
    umi_adpt_trimmed_path="${outdir_preprocess}/cutadapt/UMI_trimmed_custom"
    fastp_dir="${outdir_preprocess}/cutadapt/fastp_filter_after_UMI_trim"
    fastp_fastqc=$fastp_dir/fastqc

    if [[ -d $outdir_preprocess ]]; then
        mkdir $outdir_preprocess
    fi

    if [[ -d  ${outdir_fastqc1} ]]; then
        mkdir $outdir_fastqc1
    fi

    if [[ -d  ${outdir_cutadapt} ]]; then
        mkdir $outdir_cutadapt
    fi

    if [[ -d  ${umi_adpt_trimmed_path} ]]; then
        mkdir $umi_adpt_trimmed_path
    fi

    if [[ -d  ${fastp_dir} ]]; then
        mkdir $fastp_dir
    fi

    if [[ -d  ${fastp_fastqc} ]]; then
        mkdir $fastp_fastqc
    fi

    
    bash preprocessing_cutadapt_steps.sh \
    $indir \
    $outdir_fastqc1 \
    $outdir_cutadapt \
    $umi_adpt_trimmed_path \
    $fastp_dir \
    $fastp_fastqc > "${report_dir}/preprocessing_cutadapt_steps.out" 2>&1

    python preprocessing/cutadapt_output_parsing.py \
    "${report_dir}/preprocessing_cutadapt_steps.out" \
    "${report_dir}/cutadapt_summary.csv"

    indir=$fastp_dir

elif [[ "$cutadapt" == true ]]; then
    fastqc_dir="${outdir}/fastqc"
    fastp_dir="${outdir}/fastp"
    fastp_dir_single_samps="${outdir}/fastp/fastp_single_samples"

    bash /home/ckalk/scripts/SplitORFs/Riboseq/preprocess_data/preprocessing_steps_general.sh \
    ${indir} \
    ${fastqc_dir} \
    ${fastp_dir_single_samps}

    indir=${fastp_dir_single_samps}
fi

################################################################################
# TRANSCRIPTOMIC ALIGNMENT                                                     #
################################################################################

if [[ $transcriptomic == true && $soft_clip == true ]]; then
    bowtie2_out_dir="${outdir}/alignment_concat_transcriptome_Ignolia"
    
    source /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/bowtie2_align_k1_only_R1.sh \
     ${bowtie2_base_name} \
     no_index \
     ${fastp_dir} \
     ${bowtie2_out_dir} \
     concat_transcriptome \
     out_reports_of_runs/transcriptomic_mapping_k1_R1_norc.out \
     > out_reports_of_runs/transcriptomic_mapping_k1_R1_norc.out 2>&1


    # count soft clipping in transcriptomic alignments
    python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/analyze_soft_clipping.py ${bowtie2_out_dir}
fi


################################################################################
# TRANSCRIPTOMIC DEDUPLICATION                                                 #
################################################################################
if [[ $transcriptomic == true && $dedup == true ]]; then

    # deduplicate UMIs
    bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/deduplication/deduplicate_umi_tools.sh \
     ${UMI_indir_transcriptomic}/filtered/q10 \
     ${UMI_indir_transcriptomic}/filtered/q10/dedup \
     transcriptomic


    bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/deduplication/deduplicate_umi_tools.sh \
     ${UMI_indir_transcriptomic}/filtered \
     $umi_dedup_outdir_transcriptomic/ \
     transcriptomic


    # count soft clipping in transcriptomic alignments
    python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/analyze_soft_clipping.py $umi_dedup_outdir_transcriptomic
fi


################################################################################
# Ribowaltz                                                                    #
################################################################################

# if [ ! -d ${UMI_indir_transcriptomic}/Ribowaltz ]; then
#     mkdir ${UMI_indir_transcriptomic}/Ribowaltz
# fi

# Rscript Ribowaltz/RiboWaltz_Michi_Vlado_1_Ignolia_single_samples.R 

################################################################################
# GENOMIC ALIGNMENT                                                            #
################################################################################
star_index="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/index"
if [[ $genomic == true && $dedup == true && $soft_clip == true && $paired_reads == true ]]; then

    genome_align_dir="${outdir}/alignment_genome"
    output_star="${genome_align_dir}/STAR/only_R1"
    bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/genome_alignment_star.sh -o ${output_star} -f ${indir} -s ${star_index} \
    -a $gtf -g $genome_fasta -i -m Extend5pOfRead1 -e only_R1_


    python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
        ${output_star} \
        STAR_align_Ribo_genome.csv

elif [[ $genomic == true && $dedup == true ]]; then
    genome_align_dir="${outdir}/alignment_genome"
    output_star="${genome_align_dir}/STAR"
    echo $gtf $indir $genome_fasta $output_star $star_index
    bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/genome_alignment_star.sh -a $gtf -e "Ens_110_" -f ${indir} \
    -g $genome_fasta -m EndToEnd -o ${output_star} -s ${star_index} # -i

    python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
    ${output_star} \
    STAR_align_Ribo_genome.csv
     
fi


if [[ $genomic == true && $dedup == true ]]; then
    umi_dedup_outdir="${output_star}/deduplicated"

    # deduplicate UMIs
    source /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/deduplication/deduplicate_umi_tools.sh \
     $output_star \
     $umi_dedup_outdir

    # filter out secondary and suppl alignments
    FILES=("${umi_dedup_outdir}"/*_dedup.bam)

    for BAM in "${FILES[@]}"
    do
        samtools view -F 256 -F 2048 -q 10 -b ${BAM} > \
         "${umi_dedup_outdir}"/$(basename $BAM .bam)_filtered.bam

         samtools index "${umi_dedup_outdir}"/$(basename $BAM .bam)_filtered.bam

         samtools idxstats "${umi_dedup_outdir}"/$(basename $BAM .bam)_filtered.bam > \
        "${umi_dedup_outdir}"/$(basename $BAM .bam)_filtered_idxstats.out

        samtools stats "${umi_dedup_outdir}"/$(basename $BAM .bam)_filtered.bam > \
        "${umi_dedup_outdir}"/$(basename $BAM .bam)_filtered_stats.out

        samtools flagstat "${umi_dedup_outdir}"/$(basename $BAM .bam)_filtered.bam > \
        "${umi_dedup_outdir}"/$(basename $BAM .bam)_filtered_flagstat.out

    done


    # # analyze soft clipping of genomic deduplciated reads
    python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/analyze_soft_clipping.py $umi_dedup_outdir


    run FeatureCounts to get mapping percentages
    Rscript /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/analyze_mappings/genome_aligned_reads_biotype_counting.R

fi



# if [[ $genomic == true ]]; then
#     genome_align_dir="${outdir}/alignment_genome"
#     output_star="${genome_align_dir}/STAR"

#     if [ ! -d $genome_align_dir ]; then
#         mkdir $genome_align_dir
#     fi

#     if [ ! -d $output_star ]; then
#         mkdir $output_star
#         bash /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/STAR_Align_genomic_23_09_25.sh \
#         -i 50 "$output_star"/index \
#         $genome_fasta \
#         $gtf
#     fi


#     echo "Starting alignment against genome"

#     shopt -s nullglob
#     files=("${indir}"/*_fastp.fastq)

#     if [[ ${#files[@]} -eq 0 ]]; then
#         files=("${indir}"/*.fastq)
#     fi

#     for i in "${files[@]}"; do
#         if [[ $i == *fastp* ]]; then
#             sample_name=$(basename "$i" _fastp.fastq)
#         else
#             sample_name=$(basename "$i" .fastq)
#         fi
#         echo $i

#         bash /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/STAR_Align_genomic_23_09_25.sh \
#         -a 16 "$output_star"/index $i \
#         "$output_star"/${sample_name} \
#         EndToEnd

#         rm "$output_star"/${sample_name}*Aligned.sortedByCoord.out.bam
#         rm "$output_star"/${sample_name}*_filtered.bam
#         rm "$output_star"/*${sample_name}.bed

#         echo "===================       Sample $sample_name mapped"

#     done

#     # # deduplicate UMIs
#     # source /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/deduplication/deduplicate_umi_tools.sh \
#     #  $output_star \
#     #  $umi_dedup_outdir
# fi



if [[ $genomic == true && $riboseqc == true ]]; then
    twobit_file="/projects/splitorfs/work/reference_files/Homo_sapiens.Ensembl110.2bit"
    gtf_file="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.no.comment.gtf"
    riboseqc_outdir="${outdir}/ORFquant"

    if [[ ! -d $riboseqc_outdir ]]; then
        mkdir $riboseqc_outdir
    fi

    if [[ $dedup == true ]]; then    
        indir=${umi_dedup_outdir}
    else
        indir=${output_star}
    fi

    conda activate riboseq_qc 
    Rscript /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/ORFquant/RiboseQC.R \
        $twobit_file \
        $gtf_file \
        $riboseqc_outdir \
        $genome_fasta \
        $indir
fi



################################################################################
# THIS SHOULD NOT BE NECESSARY AS THE SO PIPELINE IS CHANGED to accomodate also 
# for dup reads, so it would be Riboseq pipeline, then SO pipeline:
# would maybe only do transcriptomic mapping and Ribowaltz for all
# but for the non-dedup samples can use the SO pipeline scripts for the SO analysis...
################################################################################
# SO ANALYSIS OF GENOMIC ALIGNMENTS                                            #
################################################################################
# bash Riboseq_SO_empirical_intersection/wrapper_empirical_intersection_random_resample.sh