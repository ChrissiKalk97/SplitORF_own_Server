#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

################################################################################
# RUN FASTQ                                                                    #
################################################################################
INDIR="/projects/splitorfs/work/own_data/Riboseq/Michi_Vlado_round_1"
OUTDIR_FASTQC1="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/fastqc_unprocessed"

# ################################################################################
# # RUN FASTQC                                                                   #
# ################################################################################
# source preprocessing/fastqc_multiqc.sh ${INDIR} ${OUTDIR_FASTQC1} fastqc_unprocessed_multiqc


################################################################################
# RUN CUTADAPT TO TRIM ADAPTERS                                                #
################################################################################
OUTDIR_CUTADAPT="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt"

# for FQ in "${INDIR}"/*R1.fastq.gz
# do
# SAMPLE=$(basename "$FQ" .R1.fastq.gz)
# cutadapt \
#     -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#     -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     -o ${OUTDIR_CUTADAPT}/${SAMPLE}.R1.cutadapt.fastq.gz \
#     -p ${OUTDIR_CUTADAPT}/${SAMPLE}.R2.cutadapt.fastq.gz \
#     --cores 16 \
#     --minimum-length 20 \
#     --max-n 0.1 \
#     --max-expected-errors 1 \
#     ${FQ} \
#     ${INDIR}/${SAMPLE}.R2.fastq.gz
# done

# Parameters:
# trim the adapter no matter where they are located and if they are parital 
# only trim at 3' end
# Require min length of 20 (inlcuding the UMI, could even be stricter)
# there is no quality filter, or possibiltiy to correct reads: plus for fastp
# --max-expected-errors could be used for quality filtering, expected errors from Phred Scores
# --max-n 0.1: filter out redas with more than 10% of N's



# ################################################################################
# # RUN FASTQC AFTER CUTADAPT                                                    #
# ################################################################################
# source preprocessing/fastqc_multiqc.sh ${OUTDIR_CUTADAPT} ${OUTDIR_CUTADAPT} cutadapt_multiqc


################################################################################
# EXTRACT UMIS CUSTOM                                                          #
################################################################################
# python preprocessing/extract_umi/extract_compare_umis.py 
# this is not optimal, think I prefer to to have the output in the same .out as the rest
# > preprocessing/extract_umi/extract_compare_umis_cutadapt.out 2>&1

Umi_adpt_trimmed_path="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/UMI_trimmed_custom"

# source preprocessing/fastqc_multiqc.sh ${Umi_adpt_trimmed_path} ${Umi_adpt_trimmed_path} UMI_trimmed_multiqc

# cd ${Umi_adpt_trimmed_path}
# rm *.fastq
# gunzip *.gz



################################################################################
# FILTERING ETC FASTP                                                          #
################################################################################
fastpOut="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim"
# source preprocessing/filter_fastp.sh ${Umi_adpt_trimmed_path} ${fastpOut} cutadapt_umi_fastp
# multiqc --force --filename ${fastpOut}/fastp_filter_after_umi_trim_multiqc ${fastpOut}





################################################################################
# Bowtie2 to align to tRNA sequences as tRNAs not present in Ensembl annotation#
################################################################################
# contamination_path="/projects/splitorfs/work/Riboseq/data/contamination"
# Bowtie_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/tRNA_aligned_bowtie2"
# source bowtie2_align.sh ${Bowtie_out_dir}/tRNA_index ${contamination_path}/hg38-tRNAs.fa ${fastpOut} ${Bowtie_out_dir} tRNA





################################################################################
# Alignment against the genome using STAR                                      #
################################################################################
OutputSTAR="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR"
source alignments/genome_alignment_star.sh ${OutputSTAR} ${fastpOut}

python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/analyze_mappings/analyze_STAR_alignments.py \
    /projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR \
    STAR_align_Ribo_genome.csv


