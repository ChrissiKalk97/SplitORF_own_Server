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
# GET CONTAMINATION SEQS                                                       #
################################################################################
# source preprocessing/contamination_filtering/prepare_Ignolia_cont_refs.sh


################################################################################
# BOWTIE ALIGN AGAINST rRNA                                                    #
################################################################################


Bowtie_base_name="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/rRNA"
Bowtie_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/contaminant_aligned"

# bowtie-build ${Bowtie_base_name}.fasta ${Bowtie_base_name}

# for FQ in "${Umi_adpt_trimmed_path}"/*R1*.fastq
# do
#     sample=$(basename "$FQ" R1.UMI_adapter_trimmed.fastq)
#     echo "${sample}"
#     FQ2="${Umi_adpt_trimmed_path}"/"${sample}"R2.UMI_adapter_trimmed.fastq
#     bowtie \
#     --best \
#     -a \
#     -y \
#     -v 3 \
#     -p 32 \
#     --allow-contain \
#     -S "${Bowtie_out_dir}"/"${sample}"rRNA_aligned.sam \
#     --un "${Bowtie_out_dir}"/"${sample}"rRNA_unaligned.fastq \
#     -x ${Bowtie_base_name} \
#     -1 ${FQ} \
#     -2 ${FQ2}
# done
# -l 22 \

# Bowtie_base_name_redownload="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA_ref_NCBI_Ens"

# bowtie-build ${Bowtie_base_name_redownload}.fasta ${Bowtie_base_name_redownload}

# for FQ in "${Umi_adpt_trimmed_path}"/*R1*.fastq
# do
#     sample=$(basename "$FQ" R1.UMI_adapter_trimmed.fastq)
#     FQ2="${Umi_adpt_trimmed_path}"/"${sample}"R2.UMI_adapter_trimmed.fastq
#     bowtie \
#     -l 22 \
#     -v 3 \
#     -p 8 \
#     --allow-contain \
#     -S "${Bowtie_out_dir}"/"${sample}"_rRNA_ref_NCBI_Ens_aligned.sam \
#     -x ${Bowtie_base_name_redownload} \
#     -1 ${FQ} \
#     -2 ${FQ2}
# done


# Bowtie2_base_name_redownload="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/bowtie2/rRNA_ref_NCBI_Ens"
# bowtie2-build ${Bowtie_base_name_redownload}.fasta ${Bowtie2_base_name_redownload} --threads 32

# for FQ in "${Umi_adpt_trimmed_path}"/*R1*.fastq
# do
#     sample=$(basename "$FQ" R1.UMI_adapter_trimmed.fastq)
#     FQ2="${Umi_adpt_trimmed_path}"/"${sample}"R2.UMI_adapter_trimmed.fastq
#     bowtie2 \
#     -p 32 \
#     -L 22 \
#     --very-sensitive-local \
#     -N 1 \
#     --no-discordant \
#     -q \
#     -S "${Bowtie_out_dir}"/"${sample}"_bowtie2_rRNA_ref_NCBI_Ens_aligned.sam \
#     -x ${Bowtie2_base_name_redownload} \
#     -1 ${FQ} \
#     -2 ${FQ2}
# done

#-N sets the number of mismatches allowed per seed, can only be 0 or 1
# not sure whether I shoudl enable this...


################################################################################
# BOWTIE ALIGN AGAINST mtrRNA                                                  #
################################################################################
# Bowtie_base_name_mtRNA="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/mtrRNArRNA45S5S"

# for FQ in "${Bowtie_out_dir}"/*rRNA_unaligned_1.fastq
# do
#     sample=$(basename "$FQ" rRNA_unaligned_1.fastq)
#     echo "${sample}"
#     FQ2="${Bowtie_out_dir}"/"${sample}"rRNA_unaligned_2.fastq
#     bowtie \
#     -a \
#     -y \
#     --best \
#     -v 3 \
#     -p 32 \
#     --allow-contain \
#     -S "${Bowtie_out_dir}"/"${sample}"_mtrRNA_aligned.sam \
#     --un "${Bowtie_out_dir}"/"${sample}"rRNA_mtrRNA_unaligned.fastq \
#     -x ${Bowtie_base_name} \
#     -1 ${FQ} \
#     -2 ${FQ2}
# done
#-l 22 \

# No more alignments found here, so not necessary to align to





################################################################################
# Bowtie2 to align to tRNA sequences as tRNAs not present in Ensembl annotation#
################################################################################
contamination_path="/projects/splitorfs/work/Riboseq/data/contamination"
Bowtie_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/tRNA_aligned_bowtie2"
source bowtie2_align.sh ${Bowtie_out_dir}/tRNA_index ${contamination_path}/hg38-tRNAs.fa ${fastpOut} ${Bowtie_out_dir} tRNA





################################################################################
# Alignment against the genome using STAR                                      #
################################################################################
OutputSTAR="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR"
source genome_alignment_star.sh ${OutputSTAR} ${Bowtie_out_dir}


