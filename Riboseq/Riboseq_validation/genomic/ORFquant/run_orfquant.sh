#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate riboseq_qc


twobit_file="/projects/splitorfs/work/reference_files/Homo_sapiens.Ensembl110.2bit"
gtf_file="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.no.comment.gtf"
outdir="/projects/splitorfs/work/Riboseq/Output/ORFquant"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
bam_path="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/NMD_genome"
for_orfquant_path="/projects/splitorfs/work/Riboseq/Output/ORFquant/RiboseQC/prepared_bam_files"

bam_path_new_huvec="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated/empirical_Riboseq_validation/NMD_genome"

export TMPDIR=/scratch/tmp/$USER



# Rscript ORFquant_prepare_annotation.R \
# $twobit_file \
# $gtf_file \
# $outdir \
# $genome_fasta \
# $bam_path

#  Rscript RiboseQC.R \
# $twobit_file \
# $gtf_file \
# $outdir \
# $genome_fasta \
# $bam_path


 Rscript RiboseQC.R \
$twobit_file \
$gtf_file \
$outdir \
$genome_fasta \
$bam_path_new_huvec

# Rscript ORFquant.R \
# $twobit_file \
# $gtf_file \
# $outdir \
# $genome_fasta \
# $for_orfquant_path

# bash intersect_ORFs_with_URs.sh \
#  /projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_NMD_16_12_24.bed \
#  /projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_RI_16_12_24.bed \
#  $for_orfquant_path
