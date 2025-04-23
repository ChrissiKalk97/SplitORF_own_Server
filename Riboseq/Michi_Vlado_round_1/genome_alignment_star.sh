#!/bin/bash

outputSTAR=$1
in_path=$2


STAR --runThreadN 50 --runMode genomeGenerate --genomeDir "$outputSTAR"/index --genomeFastaFiles /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa\
 --sjdbGTFfile /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf --sjdbOverhang 49
# --sjdbOverhang 49: maxreadlength - 1: this should actually be 37 for the leukemia samples, but for the others 30 should be fine

StarIndex="$outputSTAR"/index


for FQ in "${in_path}"/*R1*.fastq
do
    sample=$(basename "$FQ")       # remove path
    sample=${sample%%R1*}          # remove R1 and everything after
    echo ${sample}
    FQ2=${FQ/R1/R2} # substitute R1 with R2

    STAR\
    --runThreadN 16\
    --alignEndsType Extend5pOfReads12\
    --outSAMstrandField intronMotif\
    --alignIntronMin 20\
    --alignIntronMax 1000000\
    --genomeDir $StarIndex\
    --readFilesIn ${FQ} ${FQ2}\
    --twopassMode Basic\
    --seedSearchStartLmax 20\
    --outSAMattributes All\
    --alignMatesGapMax 20\
    --peOverlapNbasesMin 25\
    --outSAMtype BAM SortedByCoordinate\
    --outFileNamePrefix ${outputSTAR/}${sample}
done

# --alignMatesGapMax 20\# maximal genomic distance between mates, would like to set this to a small values as RPFs should fully overlap

# --seedSearchStartLmax 20\ reads will be split in pieces of 20 bp for MMP

# --peOverlapNbasesMin Minimum overlap of mates should be at least 30 bp, could also be set to 25 as this is the minimum
# minimum number of overlap bases to trigger mates merging and realignment

# --peOverlapMMp proportion of mismatches for overlapping sequence of the mates

# --alignEndsType Extend5pOfReads12: allow soft-clipping on 3' ends of both reads