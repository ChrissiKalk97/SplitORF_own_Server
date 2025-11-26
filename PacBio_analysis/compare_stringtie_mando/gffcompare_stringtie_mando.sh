#!/bin/bash

#----- This script runs gffcompare on 2 assemblies to merge and compare them ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio

reference_gtf=$1
mando_rescued_gtf=$2
stringtie_gtf=$3
outdir=$4
prefix=$5
genome_fasta=$6

if [[ ! -d "$outdir" ]]; then
    mkdir $outdir
fi

# gffcompare \
#     -r $reference_gtf \
#     -o $outdir/$prefix \
#     -R \
#     -K \
#     -V \
#     $mando_rescued_gtf $stringtie_gtf


#################################################################################
# ------------------ RUN SPLIT-ORFs PIPELINE                 ------------------ #
#################################################################################
bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/SplitORF_scripts/run_splitorf_pipeline_on_assembly.sh \
 $outdir/$prefix.combined.gtf \
 /home/ckalk/tools/SplitORF_pipeline \
 $genome_fasta
    