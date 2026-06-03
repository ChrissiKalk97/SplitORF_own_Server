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

python ${sqanti_path}/sqanti3_rescue.py rules \
--filter_class $classification_file \
--mode full \
--refFasta ${genome_fasta} \
--output ${output_name} \
--rescue_isoforms ${isoforms_qc} \
--rescue_gtf ${gtf_filtered} \
--refGTF ${reference_gtf} \
--refClassif ${ref_classification} \
-d ${sqanti_rescue_outdir} \
-j ${filter} \
-c 16