#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

assembly_path="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion"
huvec_assembly=$assembly_path/HUVEC/HUVEC_mando_gene_id.gtf
cm_assembly=$assembly_path/CM/CM_mando_gene_id.gtf
huvec_fasta="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/HUVEC/HUVEC_mando_gene_id_correct.fasta"
cm_fasta="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/CM/CM_mando_gene_id_correct.fasta"
genome_fasta_file="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

isoseq_reads_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"


if [[ ! -d "$assembly_path/HUVEC/pbmm2_align" ]]; then
    mkdir "$assembly_path/HUVEC/pbmm2_align"
fi

if [[ ! -d "$assembly_path/CM/pbmm2_align" ]]; then
    mkdir "$assembly_path/CM/pbmm2_align"
fi

if [[ ! -d "$assembly_path/HUVEC/pbmm2_align/HUVEC_index" ]]; then
    mkdir "$assembly_path/HUVEC/pbmm2_align/HUVEC_index"
fi

if [[ ! -d "$assembly_path/CM/pbmm2_align/CM_index" ]]; then
    mkdir "$assembly_path/CM/pbmm2_align/CM_index"
fi

if [[ ! -d "$assembly_path/CM/pbmm2_align/genome" ]]; then
    mkdir "$assembly_path/CM/pbmm2_align/genome"
fi

if [[ ! -d "$assembly_path/HUVEC/pbmm2_align/genome" ]]; then
    mkdir "$assembly_path/HUVEC/pbmm2_align/genome"
fi

if [[ ! -d "$assembly_path/HUVEC/pbmm2_align/genome_index" ]]; then
    mkdir "$assembly_path/HUVEC/pbmm2_align/genome_index"
fi


#################################################################################
# ------------------ ALIGN TO GENOME                         ------------------ #
#################################################################################

# pbmm2 index ${genome_fasta_file} $assembly_path/HUVEC/pbmm2_align/genome_index/genome_index.mmi --preset ISOSEQ

# for bam in "${isoseq_reads_dir}"/HUVEC*bam; do
#     sample=$(basename $bam)
#     sample="${sample%_merged_lima_refined.bam}"
#     pbmm2 align $assembly_path/HUVEC/pbmm2_align/genome_index/genome_index.mmi \
#     $bam \
#     $assembly_path/HUVEC/pbmm2_align/genome/${sample}_pbmm2_aligned_genome.bam \
#     --preset ISOSEQ
#  done

#  for bam in "${isoseq_reads_dir}"/CM*bam; do
#     sample=$(basename $bam)
#     sample="${sample%_merged_lima_refined.bam}"
#     pbmm2 align $assembly_path/HUVEC/pbmm2_align/genome_index/genome_index.mmi \
#     $bam \
#     $assembly_path/CM/pbmm2_align/genome/${sample}_pbmm2_aligned_genome.bam \
#     --preset ISOSEQ
#  done


# for bam in $assembly_path/HUVEC/pbmm2_align/genome/*.bam; do
#     samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
#     samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
#     samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_idxstats.out
#      samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_flagstat.out
# done 

for bam in $assembly_path/CM/pbmm2_align/genome/*.bam; do
    samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
    samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
    samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
     $(dirname $bam)/$(basename $bam .bam)_idxstats.out
     samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
     $(dirname $bam)/$(basename $bam .bam)_flagstat.out
done 

#  Rscript featureCounts_L_mode.R \
#  $assembly_path/HUVEC/pbmm2_align/genome \
#  $assembly_path/HUVEC/HUVEC_mando_gene_id.gtf

#################################################################################
# ------------------ ALIGN TO TRANSCRIPTOME                  ------------------ #
#################################################################################
# pbmm2 index ${huvec_fasta} $assembly_path/HUVEC/pbmm2_align/HUVEC_index/HUVEC_index.mmi --preset ISOSEQ
# pbmm2 index ${cm_fasta} $assembly_path/CM/pbmm2_align/CM_index/CM_index.mmi --preset ISOSEQ
# for bam in "${isoseq_reads_dir}"/HUVEC*bam; do
#     sample=$(basename $bam)
#     sample="${sample%_merged_lima_refined.bam}"
#     pbmm2 align $assembly_path/HUVEC/pbmm2_align/HUVEC_index/HUVEC_index.mmi \
#     $bam \
#     $assembly_path/HUVEC/pbmm2_align/${sample}_pbmm2_aligned.bam \
#     --preset ISOSEQ \
#     --log-level INFO
#  done

#  for bam in "${isoseq_reads_dir}"/CM*bam; do
#     sample=$(basename $bam)
#     sample="${sample%_merged_lima_refined.bam}"
#     pbmm2 align $assembly_path/CM/pbmm2_align/CM_index/CM_index.mmi \
#     $bam \
#     --preset ISOSEQ \
#     $assembly_path/CM/pbmm2_align/${sample}_pbmm2_aligned.bam \
#     --log-level INFO
#  done


# for bam in $assembly_path/HUVEC/pbmm2_align/*.bam; do
#     samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
#     samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
#     samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_idxstats.out
#      samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_flagstat.out
# done 



# for bam in $assembly_path/CM/pbmm2_align/*.bam; do
#     samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
#     samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
#     samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_idxstats.out
#      samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_flagstat.out
# done 

# python idxstats_analysis.py \
#  $assembly_path/HUVEC/pbmm2_align

#  python idxstats_analysis.py \
#  $assembly_path/CM/pbmm2_align

# python plot_venn_isoforms_by_samples.py \
#  $assembly_path/HUVEC/pbmm2_align/HUVEC_idx_summary.csv

#  python plot_venn_isoforms_by_samples.py \
#  $assembly_path/CM/pbmm2_align/CM_idx_summary.csv