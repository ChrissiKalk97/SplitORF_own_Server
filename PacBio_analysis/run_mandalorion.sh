#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

# ensembl_filtered_gtf=$1
# genome_fasta=$2
# consensus_reads_fofn=$3


ensembl_filtered_gtf="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
consensus_reads_fofn_HUVEC="pacbio_consensus_HUVEC.fofn"
consensus_reads_fofn_CM="./pacbio_consensus_CM.fofn"
out_path="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion"
bam_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"


if [[ ! -d "$out_path" ]]; then
    mkdir $out_path
fi

if [[ ! -d "$out_path/HUVEC" ]]; then
    mkdir $out_path/HUVEC
fi

if [[ ! -d "$out_path/CM" ]]; then
    mkdir $out_path/CM
fi

if [[ ! -d "$bam_dir/fastq" ]]; then
    mkdir $bam_dir/fastq
fi




shopt -s nullglob
bam_files=("${bam_dir}"/*bam)


#################################################################################
# ------------------ GET FASTQ FILES                         ------------------ #
#################################################################################
# for bam in "${bam_files[@]}"; 
# do
#     sample=$(basename $bam .bam)
#     bam2fastq -u -o $bam_dir/fastq/$sample $bam
# done

#################################################################################
# ------------------ RUN MANDOLORION TO CREATE ASSEMBLY      ------------------ #
#################################################################################
# python3 /home/ckalk/scripts/SplitORFs/PacBio_analysis/tools/Mandalorion/Mando.py -p $out_path/HUVEC -t 32 -g $ensembl_filtered_gtf -G $genome_fasta -f /projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_50NMD_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_5NMD_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_DHYPO_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_DMSO_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_DNOR_merged_lima_refined.fastq

# python3 /home/ckalk/scripts/SplitORFs/PacBio_analysis/tools/Mandalorion/Mando.py -p $out_path/CM -t 32 -g $ensembl_filtered_gtf -G $genome_fasta -f /projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/CM_5NMD_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/CM_DHYPO_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/CM_DNOR_merged_lima_refined.fastq

#################################################################################
# ------------RENAME ASSEMBLY TO GET GENE ID AND TID ENSEMBL ------------------ #
#################################################################################
# python mandalorion/rename_gene_id_name_mando_gtf.py $out_path/HUVEC/Isoforms.filtered.clean.gtf $out_path/HUVEC/HUVEC_mando_gene_id.gtf
# python mandalorion/rename_gene_id_name_mando_gtf.py $out_path/CM/Isoforms.filtered.clean.gtf $out_path/CM/CM_mando_gene_id.gtf

# the -L option does not exist and does not make any difference
# if [[ ! -d "$out_path/HUVEC/gffcompar_L_opt" ]]; then
#     mkdir $out_path/HUVEC/gffcompar_L_opt
# fi

# gffcompare -L -o $out_path/HUVEC/gffcompar_L_opt/HUVEC_gffcompar_L_opt\
#  -r $ensembl_filtered_gtf\
#   $out_path/HUVEC/HUVEC_mando_gene_id.gtf

if [[ ! -d "$out_path/HUVEC/gffcompare_only_r" ]]; then
    mkdir $out_path/HUVEC/gffcompare_only_r
fi

gffcompare -o $out_path/HUVEC/gffcompare_only_r/HUVEC_gffcompare_only_r\
 -r $ensembl_filtered_gtf\
  $out_path/HUVEC/HUVEC_mando_gene_id.gtf


  if [[ ! -d "$out_path/CM/gffcompare_only_r" ]]; then
    mkdir $out_path/CM/gffcompare_only_r
fi

gffcompare -o $out_path/CM/gffcompare_only_r/CM_gffcompare_only_r\
 -r $ensembl_filtered_gtf\
  $out_path/CM/CM_mando_gene_id.gtf

# Rename the transcripts as Ensembl IDs where applicable
python #rename script\
 #old_gtf\
  #tmap\
    #new_gtf

