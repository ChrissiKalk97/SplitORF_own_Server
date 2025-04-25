#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq


fastpOut="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/fastp_filter_after_UMI_trim"


################################################################################
# GET CONTAMINATION SEQS                                                       #
################################################################################
# source preprocessing/contamination_filtering/prepare_Ignolia_cont_refs.sh

Bowtie2_rRNA_fasta="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/rRNA.fasta"
Bowtie2_base_name="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/bowtie2/rRNA"
Bowtie2_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/Ignolia/rRNA_aligned_bowtie2"

# source bowtie2_align.sh ${Bowtie2_base_name} ${Bowtie2_rRNA_fasta} ${fastpOut} ${Bowtie2_out_dir} rRNA
# source bowtie2_align.sh ${Bowtie2_base_name} no_index ${fastpOut} ${Bowtie2_out_dir} rRNA







################################################################################
# Bowtie2 to align to tRNA sequences as tRNAs not present in Ensembl annotation#
################################################################################
tRAX_path="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/trnadb/human-trnadb-tRNAgenome"
tRAX_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/Ignolia/tRAX_aligned_bowtie2"
# source bowtie2_align.sh ${tRAX_path} no_index ${Bowtie2_out_dir} ${tRAX_out_dir} tRNA





################################################################################
# Bowtie2 to align to ncRNAs as provided by the Ignolia paper                  #
################################################################################
nc_index_path="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/ncRNAEnsembl_Gencode"
nc_fasta_file="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/Gencode_Ensembl_l_ncRNA.fa"
nc_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/Ignolia/ncRNA_aligned_bowtie2"

# however this index is for bowtie, and bowtie2 cannot use them... need to recreate the index





################################################################################
# align against mRNA MANE transcripts                                          #
################################################################################
mRNA_MANE_fasta="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna"
mRNA_align_output="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_transcriptome"
mRNA_index_path="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_transcriptome/index"


# source bowtie2_align.sh ${mRNA_index_path} ${mRNA_MANE_fasta} ${fastpOut} ${mRNA_align_output} mRNA

source bowtie2_align.sh ${mRNA_index_path} no_index ${fastpOut} ${mRNA_align_output} mRNA












################################################################################
# BOWTIE ALIGN AGAINST rRNA                                                    #
################################################################################


# Bowtie_base_name="/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/rRNA"
# Bowtie_out_dir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/cutadapt/contaminant_aligned"

# # bowtie-build ${Bowtie_base_name}.fasta ${Bowtie_base_name}

# for FQ in "${fastpOut}"/*R1*.fastq
# do
#     sample=$(basename "$FQ")       # remove path
#     sample=${sample%%R1*}          # remove R1 and everything after
#     echo ${sample}
#     FQ2=${FQ/R1/R2} # substitute R1 with R2
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







