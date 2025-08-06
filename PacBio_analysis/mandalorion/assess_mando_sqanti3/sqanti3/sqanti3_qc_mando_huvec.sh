#!/bin/bash 
eval "$(conda shell.bash hook)"
source activate sqanti3

sqanti_path=$1
gtf_file=$2
reference_gtf=$3
genome_fasta=$4
sqanti_qc_outdir=$5
short_read_file=$6
kallisto_quant_mando_raw=$7

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/gcc/x86_64-conda-linux-gnu/15.1.0:$LD_LIBRARY_PATH"

if [ ! -d "${sqanti_qc_outdir}" ]; then
    mkdir "${sqanti_qc_outdir}"
fi

# get the txt file of R1 space R2
short_read_file="sqanti3/huvec_short_reads.txt"

> "$short_read_file"  # Clear or create the output file

for r1 in $outdir_fastp/*_merged_fastp.R1.fastp.fastq.gz; do
    r2="${r1/R1/R2}"
    if [ -f "$r2" ]; then
        echo "$r1 $r2" >> "$short_read_file"
    else
        echo "Warning: No matching R2 for $r1" >&2
    fi
done

# comma separated list of kallisto quant tsv abundance files
kallisto_quant_files=$(find $kallisto_quant_mando_raw -type f -name "*.tsv" | paste -sd,)

echo $kallisto_quant_files


python ${sqanti_path}/sqanti3_qc.py \
 --isoforms ${gtf_file} \
 --refGTF ${reference_gtf} \
 --refFasta ${genome_fasta}\
 -d ${sqanti_qc_outdir} \
 --short_reads ${short_read_file} \
 --expression $kallisto_quant_files \
 --force_id_ignore \
 --polyA_motif_list ${sqanti_path}/data/polyA_motifs/mouse_and_human.polyA_motif.txt