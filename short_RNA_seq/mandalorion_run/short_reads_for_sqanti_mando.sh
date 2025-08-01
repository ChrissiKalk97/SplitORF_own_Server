#!/bin/bash

#----- This script performs the preprocessing steps for the short RNA-seq samples ----- #
# ----- and then uses them for SQANTI3 run of the raw HUVEC MAndo assembly        ----- #
# ----- FASTQC is run to obtain quality control metrics and plots                 ----- #
# ----- kallisto for quantification as required for the SQANTI3 pipeline          ----- #

eval "$(conda shell.bash hook)"
conda activate Riboseq



################################################################################
# PATH DEFINTIONS                                                              #
################################################################################
# preprocessing directories
# raw_data_dir="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData"
merged_data_dir="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/merged"
# raw_data_fastqc_dir="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData/fastqc"
# multiQC_outname="multiqc"

# reference file directories
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
gtf_file="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/HUVEC/HUVEC_mando_gene_id.gtf"
decoys="/projects/splitorfs/work/reference_files/decoys.txt"
transcript_fasta="/projects/splitorfs/work/PacBio/merged_bam_files/mandalorion/HUVEC/HUVEC_mando_gene_id_correct.fasta"
reference_gtf="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"


mando_dir_raw="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Mandalorion_raw"
outidr_fastp="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/fastp"
kallisto_quant_mando_raw="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/Mandalorion_raw/kallisto/quant"
kallisto_index_path="${mando_dir_raw}"/kallisto/index/HUVEC
sqanti_qc_outdir="${mando_dir_raw}"/SQANTI3_QC





################################################################################
# DIRECTORY CREATION                                                           #
################################################################################
if [ ! -d "${mando_dir_raw}" ]; then
    mkdir "${mando_dir_raw}"
fi

if [ ! -d "${mando_dir_raw}"/kallisto ]; then
    mkdir "${mando_dir_raw}"/kallisto
fi

if [ ! -d "${mando_dir_raw}"/kallisto/index ]; then
    mkdir "${mando_dir_raw}"/kallisto/index
fi

if [ ! -d "${mando_dir_raw}"/kallisto/index ]; then
    mkdir "${mando_dir_raw}"/SQANTI3_QC
fi


################################################################################
# PREPROCESSING STEPS                                                          #
################################################################################
# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/fastqc_multiqc_RNA_seq.sh ${raw_data_dir} ${raw_data_fastqc_dir} ${multiQC_outname} raw

# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/merge_fastq_files.sh ${raw_data_dir} ${merged_data_dir}

# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/fastqc_multiqc_RNA_seq.sh ${merged_data_dir} ${merged_data_dir}/fastqc ${multiQC_outname} raw

# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/trim_adapters_RNA_seq.sh ${merged_data_dir} ${outidr_fastp} _fastp

# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/fastqc_multiqc_RNA_seq.sh ${outidr_fastp} ${outidr_fastp}/fastqc ${multiQC_outname} fastp


################################################################################
# KALLISTO ON RAW MANDO HUVEC ASSEMBLY                                         #
################################################################################
# bash kallisto/kallisto_index.sh \
#   ${gtf_file} \
#   ${genome_fasta} \
#   ${transcript_fasta} \
#   ${kallisto_index_path}

# bash kallisto/kallisto_quantification.sh \
#  ${kallisto_index_path}.idx \
#  ${outidr_fastp} \
#  ${kallisto_quant_mando_raw}


################################################################################
# SQANTI3 QC on HUVEC assembly                                                 #
################################################################################
# create a list of all of the fastq.gz files of the short reads
# one line per sample and this needs to be R1 space R2 file
short_read_file="sqanti3/huvec_short_reads.txt"

> "$short_read_file"  # Clear or create the output file

for r1 in "$outdir_fastp"/*_merged_fastp.R2.fastp.fastq.gz; do
    r2="${r1/_R1/_R2}"
    if [ -f "$r2" ]; then
        echo "$r1 $r2" >> "$short_read_file"
    else
        echo "Warning: No matching R2 for $r1" >&2
    fi
done
 
kallisto_quant_files=$(find $kallisto_quant_mando_raw -type f -name "*.tsv" | paste -sd,)

echo $kallisto_quant_files

conda activate sqanti3
# otherwise there is an error 
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/gcc/x86_64-conda-linux-gnu/15.1.0:$LD_LIBRARY_PATH"

bash sqanti3/sqanti3_qc_mando_huvec.sh \
 /home/ckalk/tools/sqanti3 \
 ${gtf_file} \
 ${reference_gtf} \
 ${genome_fasta} \
 ${sqanti_qc_outdir} \
 sqanti3/huvec_short_reads.txt \
 ${kallisto_quant_files}







