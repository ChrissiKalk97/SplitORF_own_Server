#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq
#This script calls the aligning scripts to align the Ribo-Seq reads to the transcriptome
# and calculates the overlap with the previously determined unique regions (main pipeline) of RI and NMD
#The following are the main directories for the pipeline outputs and for the alignment outputs
genomeFasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
genomeGTF="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
outputRibotricer="/projects/splitorfs/work/Riboseq/Output/Ribotricer/URs_as_ORFs"


if [ ! -d $outputRibotricer ];then
	mkdir $outputRibotricer
fi


# Input Riboseq fastq files
input_data="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample/RI_genome/*_sorted.bam"


# First step: prepare Ribotricer ORFs
# ribotricer prepare-orfs --gtf $genomeGTF --fasta $genomeFasta --prefix "$outputRibotricer"/Ribotricer
for i in $input_data; do
sample_name=$(basename "$i" _fastp.fastq)
ribotricer  detect-orfs \
             --bam  $i \
             --ribotricer_index /projects/splitorfs/work/Riboseq/Output/Ribotricer/URs_as_ORFs/NMD_URs_as_ORFs.tsv \
             --prefix "$outputRibotricer"/${sample_name}_detected_translation

echo "Predicted translation for sample ${sample_name}"
done

