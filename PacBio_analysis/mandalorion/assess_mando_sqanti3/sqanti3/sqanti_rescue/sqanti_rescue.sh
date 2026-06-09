#!/bin/bash 
eval "$(conda shell.bash hook)"
source activate sqanti3.6

sqanti_path=$1
gtf_filtered=$2
reference_gtf=$3
genome_fasta=$4
classification_file=$5
filter=$6
output_name=$7
ref_classification=$8
isoforms_qc=$9
sqanti_rescue_outdir=${10}

echo $sqanti_rescue_outdir

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib"
echo "LD_LIBRARY_PATH is: $LD_LIBRARY_PATH"


if [ ! -d "${sqanti_rescue_outdir}" ]; then
    mkdir "${sqanti_rescue_outdir}"
fi

python ${sqanti_path}/sqanti3_rescue.py -s rules \
--filter_class $classification_file \
--mode full \
-rf ${genome_fasta} \
--output ${output_name} \
--corrected_isoforms_fasta ${isoforms_qc} \
--filtered_isoforms_gtf ${gtf_filtered} \
-rg ${reference_gtf} \
-k ${ref_classification} \
-d ${sqanti_rescue_outdir} \
-j ${filter} \
-c 16