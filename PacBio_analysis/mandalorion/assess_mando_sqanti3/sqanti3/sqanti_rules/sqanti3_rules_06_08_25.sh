#!/bin/bash 

eval "$(conda shell.bash hook)"
source activate sqanti3

sqanti_path=$1
sqanti_qc_outname=$2
sqanti_filter_outdir=$3
filter_json=$4

if [ ! -d "${sqanti_filter_outdir}" ]; then
    mkdir "${sqanti_filter_outdir}"
fi

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/gcc/x86_64-conda-linux-gnu/15.1.0:$LD_LIBRARY_PATH"


python ${sqanti_path}/sqanti3_filter.py rules \
--sqanti_class ${sqanti_qc_outname}_classification_TPM.txt \
--filter_gtf ${sqanti_qc_outname}_corrected.gtf  \
--filter_faa ${sqanti_qc_outname}_corrected.faa \
--filter_isoforms ${sqanti_qc_outname}_corrected.fasta \
-d ${sqanti_filter_outdir} \
--json ${filter_json}



