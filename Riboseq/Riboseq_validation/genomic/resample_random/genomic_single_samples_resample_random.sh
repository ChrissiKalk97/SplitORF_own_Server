#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq


outputSTAR="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10"
nmd=${outputSTAR}"/NMD_genome"
unique_region_dir="/projects/splitorfs/work/Riboseq/data/region_input/genomic"
ri=${outputSTAR}"/RI_genome"


if [ ! -d $outputSTAR ];then
	mkdir $outputSTAR
fi

if [ ! -d $nmd ];then
        mkdir $nmd
fi
if [ ! -d $ri ];then
        mkdir $ri
fi


# The following are the paths to the riboseq reads
input_data="/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples/*.fastq"

# file with the 3'UTR background regiosn
three_primes="/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed"

# Create a Logfile for the alignments in the output directory
exec > >(tee -i $outputSTAR/AlignmentLogfile.txt)
exec 2>&1

#The following block calls the Bowtie_Align script, which creates a BOWTIE index for thetranscriptome and aligns
#Ribo-seq data against it before checking the overlap with the determined unique regions. 
#For further analysis a file with random regions from the 3' and 5' UTR
#is also created and used to determine background overlap

###everything in between should not be quoted, just to be faster#############################################
echo "Starting alignment against genome"


# STAR --runThreadN 50 --runMode genomeGenerate --genomeDir "$outputSTAR"/index --genomeFastaFiles /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa\
#  --sjdbGTFfile /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf --sjdbOverhang 34
# # # --sjdbOverhang 34: maxreadlength - 1: fastp filtering of maxlength 35
# for i in $input_data; do
# sample_name=$(basename "$i" _fastp.fastq)

# source ./STAR_Align_genomic_28_01_25_resample_random_q10.sh 16 /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/index $i\
#  $nmd/${sample_name}_NMD $unique_region_dir/Unique_DNA_regions_genomic_NMD_16_12_24.bed $three_primes

# rm $nmd/*Aligned.sortedByCoord.out.bam
# rm $nmd/${sample_name}_NMD*_filtered.bam
# rm $nmd/*NMD.bed

# echo "===================       Sample $sample_name finished"

# done


# for i in $input_data; do
# sample_name=$(basename "$i" _fastp.fastq)

# source ./STAR_Align_genomic_28_01_25_resample_random_q10.sh 16 /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/index $i\
#  $ri/${sample_name}_RI $unique_region_dir/Unique_DNA_regions_genomic_RI_16_12_24.bed $three_primes

# rm $ri/*Aligned.sortedByCoord.out.bam
# rm $ri/${sample_name}_RI*_filtered.bam
# rm $ri/*RI.bed
# echo "===================       Sample $sample_name finished"

# done



export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
export MKL_ENABLE_INSTRUCTIONS=SSE4_2

# Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos="http://cran.us.r-project.org")'

# R -e 'library(rmarkdown); rmarkdown::render(input = "RiboSeqReportGenomic_iteration.Rmd", output_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/Riboseq_report_genomic_iteration_q10_16_06_25.pdf", params=list(args = c("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10", "/projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results")))'

Rscript plotting/Upsetplot_SO_upset_by_type_ISMB_talk_01_07_25.R $nmd $ri