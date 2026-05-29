#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio


# reference_gtf="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
# ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"
# genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
# consensus_reads_fofn="pacbio_consensus_${cell_type}.fofn"
# out_path="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion_updated_parameters"
# bam_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"


# kallisto_dir="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Mandalorion_raw_updated_parameters"
# short_read_dir="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/${cell_type}_fastp"
# sqanti_script_dir="/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3"



# Define the usage function
usage() {
  echo "Usage: $0 -b <bam_dir> -c <cell_type> -f <consensus_reads_fofn> -g <genome_fasta> -o <outdir_tama>  [-h for help]"
}

# Process options with silent error mode
while getopts "b:c:e:f:g:k:l:o:p:q:r:s:h" opt; do
  case $opt in
    b)
      bam_dir="$OPTARG"
      ;;
    c)
      cell_type="$OPTARG"
      ;;
    e)
      ensembl_full_gtf="$OPTARG"
      ;;
    f)
      consensus_reads_fofn="$OPTARG"
      ;;
    g)
      genome_fasta="$OPTARG"
      ;;
    k)
      kallisto_dir="$OPTARG"
      ;;
    l)
      long_read_string="$OPTARG"
      ;;
    o)
      out_path="$OPTARG"
      ;;
    p)
      mando_params="$OPTARG"
      ;;
    q)
      sqanti_script_dir="$OPTARG"
      ;;
    r) 
      reference_gtf="$OPTARG"
      ;;
    s)
      short_read_dir="$OPTARG"
      ;;
    h)
      usage
      exit 0
      ;;
    :)
      echo "Error: Option -$OPTARG requires an argument."
      usage
      exit 1
      ;;
    \?)
      echo "Error: Invalid option -$OPTARG"
      usage
      exit 1
      ;;
  esac
done











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

    if [[ -n "$mando_params" ]]; then
        IFS=';' read -r min_ratio min_reads min_feat_count upstream_buffer downstream_buffer <<< $mando_params

            python3 /home/ckalk/scripts/SplitORFs/PacBio_analysis/tools/Mandalorion/Mando.py \
            -p $out_path/${cell_type} \
            -t 32 \
            -g $reference_gtf \
            -G $genome_fasta \
            --minimum_ratio $min_ratio \
            --minimum_reads $min_reads \
            --minimum_feature_count $min_feat_count \
            -u $upstream_buffer \
            -d $downstream_buffer \
            --Acutoff 1 \
            -f ${long_read_string}
    else
    
        python3 /home/ckalk/scripts/SplitORFs/PacBio_analysis/tools/Mandalorion/Mando.py \
            -p $out_path/${cell_type} \
            -t 32 \
            -g $reference_gtf \
            -G $genome_fasta \
            --minimum_ratio 0 \
            --minimum_reads 2 \
            --minimum_feature_count 2 \
            --Acutoff 1 \
            -f ${long_read_string}
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
if [ ! -e "$out_path/${cell_type}/${cell_type}_mando_gene_id.gtf" ]; then
    python mandalorion/rename_gene_id_name_mando_gtf.py \
        $out_path/${cell_type}/Isoforms.filtered.clean.gtf \
        $out_path/${cell_type}/${cell_type}_mando_gene_id.gtf

    gffread $out_path/${cell_type}/${cell_type}_mando_gene_id.gtf \
        -g $genome_fasta -w $out_path/${cell_type}/${cell_type}_mando_gene_id.fasta
fi




################################################################################
# ------------------ LR support/expression of isoforms       ------------------ #
################################################################################
if [[ ${cell_type} == "HUVEC" ]]; then
    python mandalorion/plot_isoform_quantification_mando.py $out_path/${cell_type} 5 50
    
elif [[ ${cell_type} == "CM" ]]; then
    python mandalorion/plot_isoform_quantification_mando.py $out_path/${cell_type} 3 50
    
fi

#################################################################################
# ------------------ fl counts for SQANTI3                   ------------------ #
#################################################################################
if [[ ! -e "${out_path}/${cell_type}/${cell_type}_fl_counts.tsv" ]];then
    conda activate Riboseq
    python ${sqanti_script_dir}/sqanti3/get_fl_count_from_mando_quant_output.py \
    $out_path/${cell_type}
fi

#################################################################################
# ------------------ Run SQANTI                              ------------------ #
#################################################################################
conda activate pacbio
bash mandalorion/assess_mando_sqanti3/run_sqanti3_on_mando_one_cell_type_20_10_25.sh \
 ${cell_type} \
 $genome_fasta \
 $reference_gtf \
 $out_path \
 ${kallisto_dir} \
 $short_read_dir \
 ${sqanti_script_dir} \
 $bam_dir






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

