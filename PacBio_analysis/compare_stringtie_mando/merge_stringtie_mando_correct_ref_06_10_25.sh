#!/bin/bash

#----- This script runs gffcompare on 2 assemblies to merge and compare them ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio

reference_gtf="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"

mando_rescued_cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf"
stringtie_cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3/CM/CM_strigntie3_assembly_renamed_filtered.gtf"
outdir_cm="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/CM"
prefix_cm="CM_mando_stringtie_combined"

mando_rescued_huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf"
stringtie_huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3/HUVEC/HUVEC_strigntie3_assembly_renamed_filtered.gtf"
outdir_huvec="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/HUVEC"
prefix_huvec="HUVEC_mando_stringtie_combined"

outdir_fastp="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025"

genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"


#################################################################################
# ------------------ MERGE WITH TAMA                          ----------------- #
#################################################################################

outdir_tama="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama"

bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/tama_steps.sh \
 -r $reference_gtf \
 -f $genome_fasta \
 -o $outdir_tama \
 -s $stringtie_huvec_gtf \
 -m $mando_rescued_huvec_gtf \
 -p $outdir_fastp \
 -c HUVEC \
 -t "/home/ckalk/tools/tama"

bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/tama_steps.sh \
 -r $reference_gtf \
 -f $genome_fasta \
 -o $outdir_tama \
 -s $stringtie_cm_gtf \
 -m $mando_rescued_cm_gtf \
 -p $outdir_fastp \
 -c CM \
 -t "/home/ckalk/tools/tama"



#################################################################################
# ------------------ RUN SPLIT-ORFs PIPELINE                 ------------------ #
#################################################################################
bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
 $outdir_tama/CM/CM_merged_tama_gene_id.gtf \
 /home/ckalk/tools/SplitORF_pipeline \
 $genome_fasta


bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
 $outdir_tama/HUVEC/HUVEC_merged_tama_gene_id.gtf \
 /home/ckalk/tools/SplitORF_pipeline \
 $genome_fasta

#################################################################################
# ------------------ RUN FIFTYNT PIPELINE                    ------------------ #
#################################################################################
bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_fiftynt_on_assembly.sh \
    $outdir_tama/CM/CM_merged_tama_gene_id.gtf \
    /home/ckalk/tools/NMD_fetaure_composition \
    $genome_fasta \
    $ensembl_full_gtf \
    CM_merged_tama_50nt.csv


bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_fiftynt_on_assembly.sh \
    $outdir_tama/HUVEC/HUVEC_merged_tama_gene_id.gtf \
    /home/ckalk/tools/NMD_fetaure_composition \
    $genome_fasta \
    $ensembl_full_gtf \
    HUVEC_merged_tama_50nt.csv

#################################################################################
# ------------------ COMPARE TO ENSEMBL FULL  ASSEMBLY       ------------------ #
#################################################################################


if [[ ! -d "$outdir_tama/HUVEC/compare_Ens_full_ref" ]]; then
    mkdir $outdir_tama/HUVEC/compare_Ens_full_ref
fi

if [[ ! -d "$outdir_tama/CM/compare_Ens_full_ref" ]]; then
    mkdir $outdir_tama/CM/compare_Ens_full_ref
fi


gffcompare -o $outdir_tama/HUVEC/compare_Ens_full_ref/HUVEC_compare_full_GTF\
 -r $ensembl_full_gtf\
  $outdir_tama/HUVEC/HUVEC_merged_tama_gene_id.gtf

mv $outdir_tama/HUVEC/HUVEC_compare_full_GTF* $outdir_tama/HUVEC/compare_Ens_full_ref


gffcompare -o $outdir_tama/CM/compare_Ens_full_ref/CM_compare_full_GTF\
 -r $ensembl_full_gtf\
  $outdir_tama/CM/CM_merged_tama_gene_id.gtf

mv $outdir_tama/CM/CM_compare_full_GTF* $outdir_tama/CM/compare_Ens_full_ref

# # which isoforms have non ejcs?
python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/get_equal_ejc_isoforms.py \
 $outdir_tama/CM/compare_Ens_full_ref/CM_compare_full_GTF.CM_merged_tama_gene_id.gtf.tmap

python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/get_equal_ejc_isoforms.py \
 $outdir_tama/HUVEC/compare_Ens_full_ref/HUVEC_compare_full_GTF.HUVEC_merged_tama_gene_id.gtf.tmap


# which isoforms are novel nmd transcripts?
 python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/count_nr_novel_nmd_transcripts.py \
 /home/ckalk/tools/NMD_fetaure_composition/Output/CM_merged_tama_50nt/CM_merged_tama_50nt.csv \
 $outdir_tama/CM/compare_Ens_full_ref/CM_merged_tama_gene_id_novel_isoforms.txt \
 --assembly_type full

 python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/count_nr_novel_nmd_transcripts.py \
 /home/ckalk/tools/NMD_fetaure_composition/Output/HUVEC_merged_tama_50nt/HUVEC_merged_tama_50nt.csv \
 $outdir_tama/HUVEC/compare_Ens_full_ref/HUVEC_merged_tama_gene_id_novel_isoforms.txt \
 --assembly_type full

#################################################################################
# This can only be run when the name of the SO output folder is known   ------- #
#################################################################################
# which isoforms are novel Split-ORF transcripts?
#  python /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/count_novel_SO_transcripts.py \
#  /home/ckalk/tools/SplitORF_pipeline/Output/run_12.09.2025-17.51.04_HUVEC_tama_merged/UniqueProteinORFPairs.txt \
#  $outdir_tama/HUVEC/compare_Ens_full_ref/HUVEC_merged_tama_gene_id_novel_isoforms.txt \
#  --assembly_type full

#  python /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/count_novel_SO_transcripts.py \
#  /home/ckalk/tools/SplitORF_pipeline/Output/run_12.09.2025-14.10.14_CM_tama_merged/UniqueProteinORFPairs.txt \
#   $outdir_tama/CM/compare_Ens_full_ref/CM_merged_tama_gene_id_novel_isoforms.txt \
#  --assembly_type full


#################################################################################
# ------------------ ALIGN TO TRANSCRIPTOME                  ------------------ #
#################################################################################
# if [[ ! -d "$outdir_tama/HUVEC/minimap2" ]]; then
#     mkdir $outdir_tama/HUVEC/minimap2
# fi

# if [[ ! -d "$outdir_tama/HUVEC/minimap2/HUVEC_index" ]]; then
#     mkdir $outdir_tama/HUVEC/minimap2/HUVEC_index
# fi

# if [[ ! -d "$outdir_tama/CM/minimap2" ]]; then
#     mkdir $outdir_tama/CM/minimap2
# fi

# if [[ ! -d "$outdir_tama/CM/minimap2/CM_index" ]]; then
#     mkdir $outdir_tama/CM/minimap2/CM_index
# fi


# map-hifi unspliced genome alignment of Pacbio Long Reads
# minimap2 -x map-hifi -d $outdir_tama/HUVEC/minimap2/HUVEC_index/HUVEC_index.mmi $outdir_tama/HUVEC/HUVEC_merged_tama_gene_id.fasta
# minimap2 -x map-hifi -d $outdir_tama/CM/minimap2/CM_index/CM_index.mmi $outdir_tama/CM/CM_merged_tama_gene_id.fasta

isoseq_reads_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq"

# minimap2 align HUVEC LRs
# for fastq in "${isoseq_reads_dir}"/HUVEC*.fastq; do
#     sample=$(basename $fastq)
#     sample="${sample%_merged_lima_refined.fastq}"
#     minimap2 -a -x map-hifi $outdir_tama/HUVEC/minimap2/HUVEC_index/HUVEC_index.mmi $fastq > $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned.sam
#     samtools view -bo $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned.bam $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned.sam
#     samtools view -F 256 -F 2048 -q 10 -b $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned.bam > $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned_filtered.bam
#     samtools sort -o  $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned_sort.bam $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned_filtered.bam
#     samtools index $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned_sort.bam
#     samtools idxstats $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned_sort.bam > \
#      $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned_idxstats.out
#     samtools flagstat $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned_sort.bam > \
#      $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned_sort_flagstat.out
# done

# minimap2 align CM LRs
# for fastq in "${isoseq_reads_dir}"/CM*.fastq; do
#     sample=$(basename $fastq)
#     sample="${sample%_merged_lima_refined.fastq}"
#     minimap2 -a -x map-hifi $outdir_tama/CM/minimap2/CM_index/CM_index.mmi $fastq > $outdir_tama/CM/minimap2/${sample}_minimap2_aligned.sam
#     samtools view -bo $outdir_tama/CM/minimap2/${sample}_minimap2_aligned.bam $outdir_tama/CM/minimap2/${sample}_minimap2_aligned.sam
#     samtools view -F 256 -F 2048 -q 10 -b $outdir_tama/CM/minimap2/${sample}_minimap2_aligned.bam > $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_filtered.bam
#     samtools sort -o $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_sort.bam $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_filtered.bam
#     samtools index $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_sort.bam
#     samtools idxstats $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_sort.bam > \
#      $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_idxstats.out
#     samtools flagstat $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_sort.bam > \
#      $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_sort_flagstat.out
# done



#################################################################################
# ------------------ ALIGN TO GENOME                         ------------------ #
#################################################################################
# if [[ ! -d "$outdir_tama/HUVEC/minimap2" ]]; then
#     mkdir $outdir_tama/HUVEC/minimap2
# fi

# if [[ ! -d "$outdir_tama/HUVEC/minimap2/HUVEC_index" ]]; then
#     mkdir $outdir_tama/HUVEC/minimap2/HUVEC_index
# fi

# if [[ ! -d "$outdir_tama/CM/minimap2" ]]; then
#     mkdir $outdir_tama/CM/minimap2
# fi

# if [[ ! -d "$outdir_tama/CM/minimap2/CM_index" ]]; then
#     mkdir $outdir_tama/CM/minimap2/CM_index
# fi



#################################################################################
# - ASSESS HOW MANY NMD TRANSCRIPTS ARE PRESENT IN EACH CONDITION        ------ #
#################################################################################
# python nmd_transcripts_per_condition.py \
#  /home/ckalk/tools/NMD_fetaure_composition/Output/HUVEC_merged_tama_50nt/HUVEC_merged_tama_50nt.csv \
#  /projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama/HUVEC/minimap2

#  python nmd_transcripts_per_condition.py \
#  /home/ckalk/tools/NMD_fetaure_composition/Output/CM_merged_tama_50nt/CM_merged_tama_50nt.csv \
#  /projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama/CM/minimap2