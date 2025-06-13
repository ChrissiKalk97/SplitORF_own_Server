#!/bin/bash

# ----- This script performs Salmon quantification of short RNA-seq library         ----- #
# ----- taking the TSL and Gencode equlaity filtered Ensembl gtf as a reference     ----- #
# ----- These results can be taken to perform DEG analysis with DeSeq2 or DTE with sleuth ----- #


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


if [ ! -d $SalmonOutDir/bootstraps ]; then
    mkdir $SalmonOutDir/bootstraps
fi

shopt -s nullglob
fq_files=("${fastq_dir}"/*R1.fastp.fastq.gz)


for FQ in "${fq_files[@]}"; 
do
    SAMPLE=$(basename "$FQ")
    SAMPLE=${SAMPLE%%.*}   
    FQ2=${file/R1/R2}
    salmon quant -i ${SalmonRefDir}/index_k31 \
    -l A \
    -1 ${FQ} \
    -2 ${FQ2} \
    --validateMappings \
    -q 32 \
    -o $SalmonOutDir/bootstraps/${SAMPLE}_bootstrap \
    --seqBias --gcBias --posBias --reduceGCMemory --numBootstraps 100
done

