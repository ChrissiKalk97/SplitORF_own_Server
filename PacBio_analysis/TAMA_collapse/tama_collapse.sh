#!/bin/bash

#----- This script runs tama merge on 2 assemblies to merge them into one GTF file ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio


# Define the usage function
usage() {
  echo "Usage: $0 -r <reference_gtf> -g <genome_fasta> -o <outdir_tama>  -c <cell_type> -t <tama_tool_path> [-h for help]"
}

# Process options with silent error mode
while getopts "b:r:f:o:c:t:h" opt; do
  case $opt in
    b)
      bam_file_dir="$OPTARG"
      ;;
    r)
      reference_gtf="$OPTARG"
      ;;
    f)
      genome_fasta="$OPTARG"
      ;;
    o)
      outdir_tama="$OPTARG"
      ;;
    c)
      cell_type="$OPTARG"
      ;;
    t)
      tama_tool_path="$OPTARG"
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




# reference_gtf="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"
# genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

# bam_file_dir="/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment/HUVEC/pbmm2_align"

# outdir_tama="/projects/splitorfs/work/PacBio/merged_bam_files/tama_collapse"

# tama_tool_path="/home/ckalk/tools/tama"

#################################################################################
# ------------------ COLLAPSE WITH TAMA                          ----------------- #
#################################################################################

cell_type="${cell_type^^}"
cell_type_small="${cell_type,,}"

outdir_fastp=${outdir_fastp}/${cell_type}_fastp

if [[ ! -d "$outdir_tama" ]]; then
    mkdir $outdir_tama
fi

# if [[ ! -d "$outdir_tama/${cell_type}" ]]; then
#     mkdir $outdir_tama/${cell_type}

    # run TAMA collapse on each condition separately
    conda activate tama
    for bam in "${bam_file_dir}"/*filtered.bam; do
        echo $bam
        sample_name=$(basename $bam _pbmm2_aligned_genome_filtered.bam)

        if [[ ! -d "$outdir_tama/${cell_type}/${sample_name}" ]]; then
            mkdir $outdir_tama/${cell_type}/${sample_name}
        fi

        # python ${tama_tool_path}/tama_go/split_files/tama_mapped_sam_splitter.py $bam 10 \
        # $outdir_tama/${cell_type}/${sample_name}/${sample_name}


        # 5 bp wobble for exon ends, 50 for TES and TSS
        # use the error calculation around the splice sites (20 bp within, 5 as the threshold, from TAMA paper)
        count=0
        max_jobs=3
        for sam in $outdir_tama/${cell_type}/${sample_name}/*.sam; do
            (
            sam_name=$(basename $sam .sam)
            python ${tama_tool_path}/tama_collapse.py \
            -s ${sam} \
            -x no_cap \
            -sj sj_priority \
            -lde 5 \
            -sjt 20 \
            -m 5 \
            -a 50 \
            -z 50 \
            -f ${genome_fasta} \
            -p ${outdir_tama}/${cell_type}/${sample_name}/${sam_name} \
            -e common_ends
            ) &

            
            ((count++))
            # Every $max_jobs, wait for the batch to finish
            if (( count % max_jobs == 0 )); then
                wait
            fi
        done
        wait
    done



    # # 
    # python ${tama_tool_path}/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
    # $outdir_tama/${cell_type}/${cell_type}_merged_tama.bed \
    # $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf

# fi
  
# #################################################################################
# # ------------------ Kallisto for SQANTI QC                   ----------------- #
# #################################################################################


# if [ ! -d "$outdir_tama"/kallisto ]; then
#     mkdir "$outdir_tama"/kallisto
# fi

# if [ ! -d "$outdir_tama"/kallisto/index ]; then
#     mkdir "$outdir_tama"/kallisto/index
# fi

# if [ ! -e "$outdir_tama"/kallisto/index/${cell_type}.idx ]; then
    

#     bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/kallisto/kallisto_index.sh \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf \
#     ${genome_fasta} \
#     "$outdir_tama"/kallisto/${cell_type}_tama_merged_assembly_transcriptome.fa \
#     "$outdir_tama"/kallisto/index/${cell_type}
# fi


# if [ ! -d "$outdir_tama/kallisto/${cell_type}_quant" ]; then
#     mkdir "$outdir_tama/kallisto/${cell_type}_quant"


#     bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/kallisto/kallisto_quantification.sh \
#     "$outdir_tama/kallisto/index/${cell_type}.idx" \
#     ${outdir_fastp} \
#     "$outdir_tama/kallisto/${cell_type}_quant"
# fi



# #################################################################################
# # ------------------              SQANTI QC                   ----------------- #
# #################################################################################
# if [ ! -d "${outdir_tama}"/SQANTI3_QC ]; then
#     mkdir "${outdir_tama}"/SQANTI3_QC
# fi

# if [ ! -d "${outdir_tama}"/SQANTI3_QC/${cell_type} ]; then
#     mkdir "${outdir_tama}"/SQANTI3_QC/${cell_type}

#     conda activate pacbio

#     bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/sqanti3/sqanti3_qc_mando_huvec.sh \
#     /home/ckalk/tools/sqanti3 \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf \
#     ${reference_gtf} \
#     ${genome_fasta} \
#     "${outdir_tama}"/SQANTI3_QC/${cell_type} \
#     /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3/sqanti3/${cell_type}_short_reads.txt \
#     "$outdir_tama/kallisto/${cell_type}_quant"
    

#     python /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/get_gene_id_tama_gtf.py \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama.gtf \
#     $outdir_tama/SQANTI3_QC/${cell_type}/isoforms_classification.txt  \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama_gene_id.gtf

#     python /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/add_source_to_tama_gtf.py \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama_gene_id.gtf \
#     $outdir_tama/${cell_type}/${cell_type}_merged_tama_trans_report.txt

# fi

