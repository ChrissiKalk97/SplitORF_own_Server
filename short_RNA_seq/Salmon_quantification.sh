#!/bin/bash

# ----- This script performs Salmon quantification of short RNA-seq library         ----- #
# ----- taking the TSL and Gencode equlaity filtered Ensembl gtf as a reference     ----- #
# ----- These results can be taken to perform DEG analysis with DeSeq2 or DTE with sleuth ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio

GenomeFasta=$1
TranscriptomeAssembly=$2
SalmonRefDir=$3
ShortReadDir=$4
SalmonOutDir=$5
Decoys=$6


if [ ! -d $SalmonRefDir ]; then
    mkdir $SalmonRefDir
fi

if [ ! -d $SalmonOutDir ]; then
    mkdir $SalmonOutDir
fi




Gentrome=${SalmonRefDir}/$(basename ${TranscriptomeAssembly} .gtf)_gentrome.fa
transcriptome_fasta=${SalmonRefDir}/$(basename ${TranscriptomeAssembly} .gtf)_transcriptome.fa
if [ ! -f "$Gentrome" ]; then
    gffread $TranscriptomeAssembly -g $GenomeFasta -w $transcriptome_fasta
    cat $transcriptome_fasta $GenomeFasta > $Gentrome
fi


if [ ! -d "$SalmonRefDir"/index_k31 ]; then
    mkdir "$SalmonRefDir"/index_k31
    salmon index -t $Gentrome -i ${SalmonRefDir}/index_k31 -k 31 -d $Decoys  
fi


if [ ! -d $SalmonOutDir/Gibbs ]; then
    mkdir $SalmonOutDir/Gibbs
fi

shopt -s nullglob
fq_files=("${ShortReadDir}"/*R1.fastp.fastq.gz)


for FQ in "${fq_files[@]}"; 
do
    echo ${FQ}
    
    Sample=$(basename "$FQ")
    Sample=${Sample%%.*}   
    FQ2=${FQ/R1/R2}
    echo ${FQ2}


    salmon quant -i ${SalmonRefDir}/index_k31 \
    -l A \
    -1 ${FQ} \
    -2 ${FQ2} \
    --validateMappings \
    --threads 32 \
    -o $SalmonOutDir/Gibbs/${Sample}_Gibbs \
    --seqBias --gcBias --posBias --reduceGCMemory --numGibbsSamples 100
done

