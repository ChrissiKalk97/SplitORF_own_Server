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

# comma separated list of kallisto quant tsv abundance files
kallisto_quant_files=$(find $kallisto_quant_mando_raw -type f -name "*.tsv" | paste -sd,)

echo $kallisto_quant_files


if [ -n "${9:-}" ]; then
  fl_counts=$8
  SR_bam=$9

  python ${sqanti_path}/sqanti3_qc.py \
 --isoforms ${gtf_file} \
 --refGTF ${reference_gtf} \
 --refFasta ${genome_fasta}\
 -d ${sqanti_qc_outdir} \
 --SR_bam ${SR_bam} \
 --short_reads ${short_read_file} \
 -t 32 \
 --expression $kallisto_quant_files \
 --force_id_ignore \
  -fl $fl_counts \
 --polyA_motif_list ${sqanti_path}/data/polyA_motifs/mouse_and_human.polyA_motif.txt

elif [ -n "${8:-}" ]; then
  fl_counts=$8

  python ${sqanti_path}/sqanti3_qc.py \
 --isoforms ${gtf_file} \
 --refGTF ${reference_gtf} \
 --refFasta ${genome_fasta}\
 -d ${sqanti_qc_outdir} \
 --short_reads ${short_read_file} \
 -t 32 \
 --expression $kallisto_quant_files \
 --force_id_ignore \
  -fl $fl_counts \
 --polyA_motif_list ${sqanti_path}/data/polyA_motifs/mouse_and_human.polyA_motif.txt


else

  python ${sqanti_path}/sqanti3_qc.py \
 --isoforms ${gtf_file} \
 --refGTF ${reference_gtf} \
 --refFasta ${genome_fasta}\
 -d ${sqanti_qc_outdir} \
 --short_reads ${short_read_file} \
 -t 32 \
 --expression $kallisto_quant_files \
 --force_id_ignore \
 --polyA_motif_list ${sqanti_path}/data/polyA_motifs/mouse_and_human.polyA_motif.txt
  
fi




