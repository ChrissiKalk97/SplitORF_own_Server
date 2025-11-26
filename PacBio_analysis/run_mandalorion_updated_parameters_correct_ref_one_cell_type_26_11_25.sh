#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

# ensembl_filtered_gtf=$1
# genome_fasta=$2
# consensus_reads_fofn=$3

cell_type=$1


ensembl_filtered_gtf="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
consensus_reads_fofn="pacbio_consensus_${cell_type}.fofn"
out_path="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters"
# with the default filter settings
out_path_filter="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion"
bam_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"


if [[ ! -d "$out_path" ]]; then
    mkdir $out_path
fi


if [[ ! -d "$bam_dir/fastq" ]]; then
    mkdir $bam_dir/fastq
fi




shopt -s nullglob
bam_files=("${bam_dir}"/*bam)


#################################################################################
# ------------------ GET FASTQ FILES                         ------------------ #
#################################################################################
for bam in "${bam_files[@]}"; 
do
    sample=$(basename $bam .bam)
    if [ ! -e "$bam_dir/fastq/$sample.fastq" ];then
        bam2fastq -u -o $bam_dir/fastq/$sample $bam
    fi
done

#################################################################################
# ------------------ RUN MANDOLORION TO CREATE ASSEMBLY      ------------------ #
#################################################################################

if [[ ! -d "$out_path/${cell_type}" ]]; then
    mkdir $out_path/${cell_type}

    if [[ ${cell_type} == "HUVEC" ]]; then
        python3 /home/ckalk/scripts/SplitORFs/PacBio_analysis/tools/Mandalorion/Mando.py \
        -p $out_path/HUVEC \
        -t 32 \
        -g $ensembl_filtered_gtf \
        -G $genome_fasta \
        --minimum_ratio 0 \
        --minimum_reads 2 \
        --minimum_feature_count 2 \
        --Acutoff 1 \
        -f /projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_50NMD_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_5NMD_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_DHYPO_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_DMSO_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/HUVEC_DNOR_merged_lima_refined.fastq

    elif [[ ${cell_type} == "CM" ]]; then
        python3 /home/ckalk/scripts/SplitORFs/PacBio_analysis/tools/Mandalorion/Mando.py \
            -p $out_path/CM \
            -t 32 \
            -g $ensembl_filtered_gtf \
            -G $genome_fasta \
            --minimum_ratio 0 \
            --minimum_reads 2 \
            --minimum_feature_count 2 \
            --Acutoff 1 \
            -f /projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/CM_5NMD_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/CM_DHYPO_merged_lima_refined.fastq,/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine/fastq/CM_DNOR_merged_lima_refined.fastq
    fi

fi

if [ ! -e "$out_path/${cell_type}/tmp/Isoforms.aligned.out.filtered.bam" ];then
    # get bam files of minimap2 alignments made with Mando
    samtools view -bo $out_path/${cell_type}/tmp/Isoforms.aligned.out.filtered.bam \
    $out_path/${cell_type}/tmp/Isoforms.aligned.out.filtered.sam 

    samtools sort $out_path/${cell_type}/tmp/Isoforms.aligned.out.filtered.bam \
    -o $out_path/${cell_type}/tmp/Isoforms.aligned.out.filtered.sorted.bam

    samtools index $out_path/${cell_type}/tmp/Isoforms.aligned.out.filtered.sorted.bam
fi


#################################################################################
# ------------RENAME ASSEMBLY TO GET GENE ID AND TID ENSEMBL ------------------ #
#################################################################################
python mandalorion/rename_gene_id_name_mando_gtf.py \
    $out_path/${cell_type}/Isoforms.filtered.clean.gtf \
    $out_path/${cell_type}/${cell_type}_mando_gene_id.gtf

gffread $out_path/${cell_type}/${cell_type}_mando_gene_id.gtf \
    -g $genome_fasta -w $out_path/${cell_type}/${cell_type}_mando_gene_id.fasta


################################################################################
# ------------------ LR support/expression of isoforms       ------------------ #
################################################################################
python mandalorion/plot_isoform_quantification_mando.py $out_path/${cell_type} 5 50
python mandalorion/plot_isoform_quantification_mando.py $out_path_filter/${cell_type} 5 50


#################################################################################
# ------------------ fl counts for SQANTI3                   ------------------ #
#################################################################################
conda activate Riboseq
python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/sqanti3/get_fl_count_from_mando_quant_output.py \
  $out_path/${cell_type}

#################################################################################
# ------------------ Run SQANTI                              ------------------ #
#################################################################################
conda activate pacbio
bash mandalorion/assess_mando_sqanti3/run_sqanti3_on_mando_one_cell_type_20_10_25.sh \
 ${cell_type} \
 $genome_fasta \
 $ensembl_filtered_gtf \
 $out_path \
 "/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Mandalorion_raw_updated_parameters" \
 "/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/${cell_type}_fastp" \
 "/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3" \
 "/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"



#################################################################################
# ------------------ RUN SPLIT-ORFs PIPELINE                 ------------------ #
#################################################################################
# bash SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
# $out_path/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf \
# /home/ckalk/tools/SplitORF_pipeline \
# $genome_fasta


# bash SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
# $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf \
# /home/ckalk/tools/SplitORF_pipeline \
# $genome_fasta

#################################################################################
# ------------------ RUN FIFTYNT PIPELINE                    ------------------ #
#################################################################################
# bash SplitORF_scripts/run_fiftynt_on_assembly.sh \
#     $out_path/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf \
#     /home/ckalk/tools/NMD_fetaure_composition \
#     $genome_fasta \
#     $ensembl_full_gtf \
#     CM_mando_rescued_50nt.csv


# bash SplitORF_scripts/run_fiftynt_on_assembly.sh \
#     $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf \
#     /home/ckalk/tools/NMD_fetaure_composition \
#     $genome_fasta \
#     $ensembl_full_gtf \
#     HUVEC_mando_rescued_50nt.csv






#################################################################################
# ------------------ COMPARE TO ENSEMBL FULL  ASSEMBLY       ------------------ #
#################################################################################


# if [[ ! -d "$out_path/SQANTI3/SQANTI3_Rescue/HUVEC/compare_Ens_full_ref" ]]; then
#     mkdir $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/compare_Ens_full_ref
# fi

# if [[ ! -d "$out_path/SQANTI3/SQANTI3_Rescue/CM/compare_Ens_full_ref" ]]; then
#     mkdir $out_path/SQANTI3/SQANTI3_Rescue/CM/compare_Ens_full_ref
# fi


# gffcompare -o $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/compare_Ens_full_ref/HUVEC_rescue_compare_full_GTF\
#  -r $ensembl_full_gtf\
#   $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_rules_filter_rescued.gtf

# mv $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/HUVEC_rescue_compare_full_GTF* $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/compare_Ens_full_ref


# gffcompare -o $out_path/SQANTI3/SQANTI3_Rescue/CM/compare_Ens_full_ref/CM_rescue_compare_full_GTF\
#  -r $ensembl_full_gtf\
#   $out_path/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_rules_filter_rescued.gtf

# mv $out_path/SQANTI3/SQANTI3_Rescue/CM/CM_rescue_compare_full_GTF* $out_path/SQANTI3/SQANTI3_Rescue/CM/compare_Ens_full_ref

# # # which isoforms have non ejcs?
# python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/get_equal_ejc_isoforms.py \
#  $out_path/SQANTI3/SQANTI3_Rescue/CM/compare_Ens_full_ref/CM_rescue_compare_full_GTF.CM_rescue_rules_filter_rescued.gtf.tmap

# python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/get_equal_ejc_isoforms.py \
#  $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/compare_Ens_full_ref/HUVEC_rescue_compare_full_GTF.HUVEC_rescue_rules_filter_rescued.gtf.tmap


# which isoforms are novel nmd transcripts?
#  python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/count_nr_novel_nmd_transcripts.py \
#  /home/ckalk/tools/NMD_fetaure_composition/Output/CM_mando_rescued_50nt/CM_mando_rescued_50nt.csv  \
#  $out_path/SQANTI3/SQANTI3_Rescue/CM/compare_Ens_full_ref/CM_rescue_rules_filter_rescued_novel_isoforms.txt \
#  --assembly_type full

#  python /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/count_nr_novel_nmd_transcripts.py \
#  /home/ckalk/tools/NMD_fetaure_composition/Output/HUVEC_mando_rescued_50nt/HUVEC_mando_rescued_50nt.csv  \
#  $out_path/SQANTI3/SQANTI3_Rescue/HUVEC/compare_Ens_full_ref/HUVEC_rescue_rules_filter_rescued_novel_isoforms.txt \
#  --assembly_type full

