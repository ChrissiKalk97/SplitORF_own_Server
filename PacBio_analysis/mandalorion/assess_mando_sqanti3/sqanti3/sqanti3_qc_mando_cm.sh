#!/bin/bash 
eval "$(conda shell.bash hook)"
source activate sqanti3

sqanti_path=$1
gtf_file=$2
reference_gtf=$3
genome_fasta=$4
sqanti_qc_outdir=$5


export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/gcc/x86_64-conda-linux-gnu/15.1.0:$LD_LIBRARY_PATH"

if [ ! -d "${sqanti_qc_outdir}"/CM ]; then
    mkdir "${sqanti_qc_outdir}"/CM
fi

python ${sqanti_path}/sqanti3_qc.py \
 --isoforms ${gtf_file} \
 --refGTF ${reference_gtf} \
 --refFasta ${genome_fasta}\
 -d ${sqanti_qc_outdir} \
 --force_id_ignore \
 --polyA_motif_list ${sqanti_path}/data/polyA_motifs/mouse_and_human.polyA_motif.txt