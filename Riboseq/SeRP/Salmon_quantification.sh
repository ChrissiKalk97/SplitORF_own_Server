#!/bin/bash

# ----- This script performs Salmon quantification of the R1's of a Riboseq library ----- #
# ----- taking the TSL and Gencode equlaity filtered Ensembl gtf as a reference     ----- #
# ----- These results can be taken to perform DEG analysis with DeSeq2              ----- #


GenomeFasta=$1
TranscriptomeAssembly=$2
SalmonRefDir=$3
ShortReadDir=$4
SalmonOutDir=$5


if [ ! -d $SalmonRefDir ]; then
    mkdir $SalmonRefDir
fi

if [ ! -d $SalmonOutDir ]; then
    mkdir $SalmonOutDir
fi




Gentrome=${SalmonRefDir}/$(basename ${TranscriptomeAssembly} .gtf)_gentrome.fa
if [ ! -f "$Gentrome" ]; then
    gffread $TranscriptomeAssembly -g $GenomeFasta -w $transcriptome_fasta
    cat $transcriptome_fasta $GenomeFasta > $Gentrome
fi


if [ ! -d "$SalmonRefDir"/index_k17 ]; then
    mkdir "$SalmonRefDir"/index_k17
    salmon index -t $Gentrome -i ${SalmonRefDir}/index_k17 -k 17 -d ${SalmonRefDir}/decoys.txt  
fi



for FQ in ${ShortReadDir}/*R1*.fastq.gz
do
    sample=$(basename "$FQ")
    sample=${sample%%.*}  
    salmon quant -i ${SalmonRefDir}/index_k17\
    -l A -r ${FQ} \
    --validateMappings\
    -q 32 \
    -o $SalmonOutDir/bootstraps/${sample}_bootstrap\
    --seqBias --gcBias --posBias --reduceGCMemory --numBootstraps 100
done

