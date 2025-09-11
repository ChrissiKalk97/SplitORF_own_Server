#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

pacbio_merged_bamdir="/projects/splitorfs/work/PacBio/merged_bam_files"

lima_outdir=${pacbio_merged_bamdir}/lima
refine_outdir=${pacbio_merged_bamdir}/isoseq/refine


echo "Raw read counts per bam"
for bam in $pacbio_merged_bamdir/*.bam; do
    echo $bam
    samtools view -c $bam
done


echo "Read counts after lima"
for bam in $lima_outdir/*.bam; do
    echo $bam
    samtools view -c $bam
done

echo "Read counts after refine"
for bam in $refine_outdir/*.bam; do
    echo $bam
    samtools view -c $bam
done