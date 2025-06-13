#!/bin/bash

fastq_dir=$1
STAR_out_dir=$2
STAR_index=$3
genome_fasta=$4
genome_annotation=$5

if [ ! -d ${STAR_out_dir} ]; then
    mkdir $STAR_out_dir
fi

if [ ! -d ${STAR_out_dir}/index ]; then
    mkdir $STAR_out_dir/index
fi


STAR \
 --runThreadN 16 \
 --runMode genomeGenerate \
 --genomeDir $STAR_index\
 --genomeFastaFiles $genome_fasta \
 --sjdbGTFfile $genome_annotation \
 --sjdbOverhang 149


shopt -s nullglob
fq_files=("${fastq_dir}"/*R1.fastp.fastq.gz)


for FQ in "${fq_files[@]}"; 
do
        SAMPLE=$(basename "$FQ")
        SAMPLE=${SAMPLE%%R1*}   
        FQ2=${file/R1/R2}
        STAR \
        --runThreadN 16 \
        --genomeDir $STAR_index\
        --readFilesIn $FQ $FQ2\
        --readFilesCommand gunzip -c \
        --alignEndsType EndToEnd \
        --alignSJoverhangMin 10 \
        --alignSJDBoverhangMin 5 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.05 \
        --outFilterMatchNminOverLread 0.7 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMstrandField intronMotif \
        --twopassMode Basic \
        --outSAMattributes NH HI AS nM NM MD jM jI XS \
        --outFileNamePrefix $STAR_out_dir/${SAMPLE}  
        
done

# wait




# --outSAMstrandField intronMotif required for Stringtie, XS tag
# --alignMatesGapMax 1000000 id the default
# --alignIntronMax 1000000 is the default
# --alignIntronMin 20 is the default
# --outFilterMatchNminOverLread 0.7 default 0.66, percentage of read that needs to align (matching bases)
# -outFilterMismatchNoverLmax 0.05, default 0.3: alignment is only output if the ratio of mismatches
# to aligned positions is smaller than this value
# --outFilterMismatchNmax 999 default, max nr of mismatches per pair, large number switches this off
# --outFilterScoreMin 1, defualt is 0, only output alignments higher than this score
# `--alignSJDBoverhangMin` default 1: minimum overahng of annotated junctions
#--alignSJoverhangMin default 8, min overahng of unannotated junctions
# --outSJfilterOverhangMin default 30 12 12 12, for unannotated jcts, minimum overhang length for splice junctions on both sides for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.