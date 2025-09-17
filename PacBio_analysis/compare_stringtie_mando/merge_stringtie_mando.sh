#!/bin/bash

#----- This script runs gffcompare on 2 assemblies to merge and compare them ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio

reference_gtf="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"
ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"

mando_rescued_cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf"
stringtie_cm_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3/CM/CM_strigntie3_assembly_renamed_filtered.gtf"
outdir_cm="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/CM"
prefix_cm="CM_mando_stringtie_combined"

mando_rescued_huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf"
stringtie_huvec_gtf="/projects/splitorfs/work/PacBio/merged_bam_files/stringtie3/HUVEC/HUVEC_strigntie3_assembly_renamed_filtered.gtf"
outdir_huvec="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/HUVEC"
prefix_huvec="HUVEC_mando_stringtie_combined"

genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"


#################################################################################
# ------------------ MERGE WITH GFFCOMPARE                   ------------------ #
#################################################################################

# bash gffcompare_stringtie_mando.sh \
#     $reference_gtf \
#     $mando_rescued_cm_gtf \
#     $stringtie_cm_gtf \
#     $outdir_cm \
#     $prefix_cm \
#     $genome_fasta


# bash gffcompare_stringtie_mando.sh \
#     $reference_gtf \
#     $mando_rescued_huvec_gtf \
#     $stringtie_huvec_gtf \
#     $outdir_huvec \
#     $prefix_huvec \
#     $genome_fasta



# script_dir="/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3"
# genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_cm/CM_mando_stringtie_combined.combined.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_cm/SQANTI_QC

#  bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_huvec/HUVEC_mando_stringtie_combined.combined.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_huvec/SQANTI_QC



#################################################################################
# ------------------ MERGE WITH STRINGTIE MERGE               ----------------- #
#################################################################################
# outdir_stringtie_merge="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/stringtie_merge"

# if [[ ! -d "$outdir_stringtie_merge" ]]; then
#     mkdir $outdir_stringtie_merge
# fi

# if [[ ! -d "$outdir_stringtie_merge/CM" ]]; then
#     mkdir $outdir_stringtie_merge/CM
# fi

# if [[ ! -d "$outdir_stringtie_merge/HUVEC" ]]; then
#     mkdir $outdir_stringtie_merge/HUVEC
# fi


# /home/ckalk/tools/stringtie-3.0.1.Linux_x86_64/stringtie \
#   --merge \
#   -o $outdir_stringtie_merge/CM/CM_mando_stringtie_merged.gtf\
#     $mando_rescued_cm_gtf \
#     $stringtie_cm_gtf

# /home/ckalk/tools/stringtie-3.0.1.Linux_x86_64/stringtie \
#   --merge \
#   -o $outdir_stringtie_merge/HUVEC/HUVEC_mando_stringtie_merged.gtf \
#     $mando_rescued_huvec_gtf \
#     $stringtie_huvec_gtf



# bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_stringtie_merge/CM/CM_mando_stringtie_merged.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_stringtie_merge/CM/SQANTI_QC

#  bash ${script_dir}/sqanti3/sqanti3_qc_mando_cm.sh \
#  /home/ckalk/tools/sqanti3 \
#  $outdir_stringtie_merge/HUVEC/HUVEC_mando_stringtie_merged.gtf \
#  ${reference_gtf} \
#  ${genome_fasta} \
#  $outdir_stringtie_merge/HUVEC/SQANTI_QC

#################################################################################
# ------------------ MERGE WITH TAMA                          ----------------- #
#################################################################################

outdir_tama="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama"

bash tama_steps.sh \
 -r $reference_gtf \
 -f $genome_fasta \
 -o $outdir_tama \
 -s $stringtie_cm_gtf \
 -m $mando_rescued_cm_gtf \
 -c CM \
 -t "/home/ckalk/tools/tama"

bash tama_steps.sh \
 -r $reference_gtf \
 -f $genome_fasta \
 -o $outdir_tama \
 -s $stringtie_huvec_gtf \
 -m $mando_rescued_huvec_gtf \
 -c HUVEC \
 -t "/home/ckalk/tools/tama"


# #################################################################################
# # ------------------ RUN SPLIT-ORFs PIPELINE                 ------------------ #
# #################################################################################
# bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
#  $outdir_tama/CM/CM_merged_tama_gene_id.gtf \
#  /home/ckalk/tools/SplitORF_pipeline \
#  $genome_fasta


# bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
#  $outdir_tama/HUVEC/HUVEC_merged_tama_gene_id.gtf \
#  /home/ckalk/tools/SplitORF_pipeline \
#  $genome_fasta

# #################################################################################
# # ------------------ RUN FIFTYNT PIPELINE                    ------------------ #
# #################################################################################
# bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_fiftynt_on_assembly.sh \
#     $outdir_tama/CM/CM_merged_tama_gene_id.gtf \
#     /home/ckalk/tools/NMD_fetaure_composition \
#     $genome_fasta \
#     $ensembl_full_gtf \
#     CM_merged_tama_50nt.csv


# bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_fiftynt_on_assembly.sh \
#     $outdir_tama/HUVEC/HUVEC_merged_tama_gene_id.gtf \
#     /home/ckalk/tools/NMD_fetaure_composition \
#     $genome_fasta \
#     $ensembl_full_gtf \
#     HUVEC_merged_tama_50nt.csv

#################################################################################
# ------------------ COMPARE TO ENSEMBL FULL  ASSEMBLY       ------------------ #
#################################################################################


# if [[ ! -d "$outdir_tama/HUVEC/compare_Ens_full_ref" ]]; then
#     mkdir $outdir_tama/HUVEC/compare_Ens_full_ref
# fi

# if [[ ! -d "$outdir_tama/CM/compare_Ens_full_ref" ]]; then
#     mkdir $outdir_tama/CM/compare_Ens_full_ref
# fi


# gffcompare -o $outdir_tama/HUVEC/compare_Ens_full_ref/HUVEC_compare_full_GTF\
#  -r $ensembl_full_gtf\
#   $outdir_tama/HUVEC/HUVEC_merged_tama_gene_id.gtf

# mv $outdir_tama/HUVEC/HUVEC_compare_full_GTF* $outdir_tama/HUVEC/compare_Ens_full_ref


# gffcompare -o $outdir_tama/CM/compare_Ens_full_ref/CM_compare_full_GTF\
#  -r $ensembl_full_gtf\
#   $outdir_tama/CM/CM_merged_tama_gene_id.gtf

# mv $outdir_tama/CM/CM_compare_full_GTF* $outdir_tama/CM/compare_Ens_full_ref

# # # which isoforms have non ejcs?
# python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/get_equal_ejc_isoforms.py \
#  $outdir_tama/CM/compare_Ens_full_ref/CM_compare_full_GTF.CM_merged_tama_gene_id.gtf.tmap

# python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/get_equal_ejc_isoforms.py \
#  $outdir_tama/HUVEC/compare_Ens_full_ref/HUVEC_compare_full_GTF.HUVEC_merged_tama_gene_id.gtf.tmap


# # which isoforms are novel nmd transcripts?
#  python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/count_nr_novel_nmd_transcripts.py \
#  /home/ckalk/tools/NMD_fetaure_composition/Output/CM_merged_tama_50nt/CM_merged_tama_50nt.csv \
#  $outdir_tama/CM/compare_Ens_full_ref/CM_merged_tama_gene_id_novel_isoforms.txt \
#  --assembly_type full

#  python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/count_nr_novel_nmd_transcripts.py \
#  /home/ckalk/tools/NMD_fetaure_composition/Output/HUVEC_merged_tama_50nt/HUVEC_merged_tama_50nt.csv \
#  $outdir_tama/HUVEC/compare_Ens_full_ref/HUVEC_merged_tama_gene_id_novel_isoforms.txt \
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
#     samtools sort -o  $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned_sort.bam $outdir_tama/HUVEC/minimap2/${sample}_minimap2_aligned.bam
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
#     samtools sort -o $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_sort.bam $outdir_tama/CM/minimap2/${sample}_minimap2_aligned.bam
#     samtools index $outdir_tama/CM/minimap2/${sample}_minimap2_aligned.bam
#     samtools idxstats $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_sort.bam > \
#      $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_idxstats.out
#     samtools flagstat $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_sort.bam > \
#      $outdir_tama/CM/minimap2/${sample}_minimap2_aligned_sort_flagstat.out
# done



