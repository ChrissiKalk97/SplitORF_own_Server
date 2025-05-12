#!/bin/bash

outputSTAR=$1
in_path=$2
StarIndex=$3

# STAR --runThreadN 50 --runMode genomeGenerate --genomeDir $Star_index --genomeFastaFiles /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa\
#  --sjdbGTFfile /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf --sjdbOverhang 49
# --sjdbOverhang 49: maxreadlength - 1: this should actually be 37 for the leukemia samples, but for the others 30 should be fine

shopt -s nullglob  # Prevents literal pattern if no match
files=(*.fatsq)

if [ ${#files[@]} -gt 0 ]; then
    cd $in_path
    gunzip *.gz
    cd -
fi

# First try matching *R1* files
shopt -s nullglob  # So non-matching globs result in empty arrays
files=("${in_path}"/*R1*.fastq)

if [ ${#files[@]} -eq 0 ]; then
    files=("${in_path}"/*.fastq.1)
fi



for FQ in "${files[@]}"
do
    sample=$(basename "$FQ")       # remove path
    if [[ ${sample} == *"R1"* ]]; then
        sample=${sample%%R1*}          # remove R1 and everything after
        FQ2=${FQ/R1/R2} # substitute R1 with R2 in the whole file path
    else
        sample=${sample%%.fastq*} 
        FQ2=${FQ/.1/.2} # substitute R1 with R2 in the whole file path
    fi
    echo ${sample}
    echo ${FQ}
    echo ${FQ2}

    STAR\
    --runThreadN 16\
    --alignEndsType Extend5pOfRead1\
    --outSAMstrandField intronMotif\
    --alignIntronMin 20\
    --alignIntronMax 1000000\
    --genomeDir $StarIndex\
    --readFilesIn ${FQ}\
    --twopassMode Basic\
    --seedSearchStartLmax 20\
    --outFilterMatchNminOverLread 0.9\
    --outSAMattributes All\
    --outSAMtype BAM SortedByCoordinate\
    --outFileNamePrefix "${outputSTAR}"/"${sample}"only_R1_

    samtools index --threads 32 "${outputSTAR}"/"${sample}"only_R1_Aligned.sortedByCoord.out.bam

    samtools idxstats "${outputSTAR}"/"${sample}"only_R1_Aligned.sortedByCoord.out.bam > \
    "${outputSTAR}"/"${sample}"only_R1_idxstats.out

    samtools stats "${outputSTAR}"/"${sample}"only_R1_Aligned.sortedByCoord.out.bam > \
    "${outputSTAR}"/"${sample}"only_R1_stats.out
done

# --alignMatesGapMax 20\# maximal genomic distance between mates, would like to set this to a small values as RPFs should fully overlap

# --seedSearchStartLmax 20\ reads will be split in pieces of 20 bp for MMP

# --peOverlapNbasesMin Minimum overlap of mates should be at least 30 bp, could also be set to 25 as this is the minimum
# minimum number of overlap bases to trigger mates merging and realignment

# --peOverlapMMp proportion of mismatches for overlapping sequence of the mates

# --alignEndsType Extend5pOfReads12: allow soft-clipping on 3' ends of both reads

# --outFilterMultimapNmax: defualt is 10, allow more in order to keep tRNA and rRNA alignments
#     --outFilterMultimapNmax 20\ makes files big and deduplication slow