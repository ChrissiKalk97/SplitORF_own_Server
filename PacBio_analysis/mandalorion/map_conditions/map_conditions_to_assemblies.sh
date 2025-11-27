#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio


out_path="/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment"
huvec_assembly=$out_path/HUVEC/HUVEC_mando_gene_id.gtf
cm_assembly=$out_path/CM/CM_mando_gene_id.gtf
huvec_fasta="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/HUVEC/HUVEC_mando_gene_id_correct.fasta"
cm_fasta="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/CM/CM_mando_gene_id_correct.fasta"
genome_fasta_file="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

isoseq_reads_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"

if [[ ! -d "$out_path" ]]; then
    mkdir "$out_path"
fi

if [[ ! -d "$out_path/HUVEC" ]]; then
    mkdir "$out_path/HUVEC"
fi

if [[ ! -d "$out_path/CM" ]]; then
    mkdir "$out_path/CM"
fi

if [[ ! -d "$out_path/HUVEC/pbmm2_align" ]]; then
    mkdir "$out_path/HUVEC/pbmm2_align"
fi

if [[ ! -d "$out_path/CM/pbmm2_align" ]]; then
    mkdir "$out_path/CM/pbmm2_align"
fi

# if [[ ! -d "$out_path/HUVEC/pbmm2_align/HUVEC_index" ]]; then
#     mkdir "$out_path/HUVEC/pbmm2_align/HUVEC_index"
# fi

# if [[ ! -d "$out_path/CM/pbmm2_align/CM_index" ]]; then
#     mkdir "$out_path/CM/pbmm2_align/CM_index"
# fi



if [[ ! -d "$out_path/HUVEC/pbmm2_align/genome_index" ]]; then
    mkdir "$out_path/HUVEC/pbmm2_align/genome_index"
fi


#################################################################################
# ------------------ ALIGN TO GENOME                         ------------------ #
#################################################################################
if [ ! -e "$out_path/HUVEC/pbmm2_align/genome_index/genome_index.mmi" ]; then 
    pbmm2 index ${genome_fasta_file} $out_path/HUVEC/pbmm2_align/genome_index/genome_index.mmi --preset ISOSEQ
fi


for bam in "${isoseq_reads_dir}"/HUVEC*bam; do
    sample=$(basename $bam)
    sample="${sample%_merged_lima_refined.bam}"
    pbmm2 align $out_path/HUVEC/pbmm2_align/genome_index/genome_index.mmi \
    $bam \
    $out_path/HUVEC/pbmm2_align/${sample}_pbmm2_aligned_genome.bam \
    --preset ISOSEQ
 done

 for bam in "${isoseq_reads_dir}"/CM*bam; do
    sample=$(basename $bam)
    sample="${sample%_merged_lima_refined.bam}"
    pbmm2 align $out_path/HUVEC/pbmm2_align/genome_index/genome_index.mmi \
    $bam \
    $out_path/CM/pbmm2_align/${sample}_pbmm2_aligned_genome.bam \
    --preset ISOSEQ
 done


for bam in $out_path/HUVEC/pbmm2_align/*.bam; do
    samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
    samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
    samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
     $(dirname $bam)/$(basename $bam .bam)_idxstats.out
     samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
     $(dirname $bam)/$(basename $bam .bam)_flagstat.out
done 

for bam in $out_path/CM/pbmm2_align/*.bam; do
    samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
    samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
    samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
     $(dirname $bam)/$(basename $bam .bam)_idxstats.out
     samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
     $(dirname $bam)/$(basename $bam .bam)_flagstat.out
done 

#  Rscript featureCounts_L_mode.R \
#  $out_path/HUVEC/pbmm2_align/genome \
#  $out_path/HUVEC/HUVEC_mando_gene_id.gtf

#################################################################################
# ------------------ ALIGN TO TRANSCRIPTOME                  ------------------ #
#################################################################################
# pbmm2 index ${huvec_fasta} $out_path/HUVEC/pbmm2_align/HUVEC_index/HUVEC_index.mmi --preset ISOSEQ
# pbmm2 index ${cm_fasta} $out_path/CM/pbmm2_align/CM_index/CM_index.mmi --preset ISOSEQ
# for bam in "${isoseq_reads_dir}"/HUVEC*bam; do
#     sample=$(basename $bam)
#     sample="${sample%_merged_lima_refined.bam}"
#     pbmm2 align $out_path/HUVEC/pbmm2_align/HUVEC_index/HUVEC_index.mmi \
#     $bam \
#     $out_path/HUVEC/pbmm2_align/${sample}_pbmm2_aligned.bam \
#     --preset ISOSEQ \
#     --log-level INFO
#  done

#  for bam in "${isoseq_reads_dir}"/CM*bam; do
#     sample=$(basename $bam)
#     sample="${sample%_merged_lima_refined.bam}"
#     pbmm2 align $out_path/CM/pbmm2_align/CM_index/CM_index.mmi \
#     $bam \
#     --preset ISOSEQ \
#     $out_path/CM/pbmm2_align/${sample}_pbmm2_aligned.bam \
#     --log-level INFO
#  done


# for bam in $out_path/HUVEC/pbmm2_align/*.bam; do
#     samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
#     samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
#     samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_idxstats.out
#      samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_flagstat.out
# done 



# for bam in $out_path/CM/pbmm2_align/*.bam; do
#     samtools sort -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
#     samtools index $(dirname $bam)/$(basename $bam .bam)_sorted.bam
#     samtools idxstats $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_idxstats.out
#      samtools flagstat $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
#      $(dirname $bam)/$(basename $bam .bam)_flagstat.out
# done 

# python idxstats_analysis.py \
#  $out_path/HUVEC/pbmm2_align

#  python idxstats_analysis.py \
#  $out_path/CM/pbmm2_align

# python plot_venn_isoforms_by_samples.py \
#  $out_path/HUVEC/pbmm2_align/HUVEC_idx_summary.csv

#  python plot_venn_isoforms_by_samples.py \
#  $out_path/CM/pbmm2_align/CM_idx_summary.csv