#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq
#This script calls the aligning scripts to align the Ribo-Seq reads to the transcriptome
# and calculates the overlap with the previously determined unique regions (main pipeline) of RI and NMD
#The following are the main directories for the pipeline outputs and for the alignment outputs
outputSTAR="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples"
nmd=${outputSTAR}"/NMD_transcriptome"
unique_region_dir="/projects/splitorfs/work/Riboseq/data/region_input/genomic"
ri=${outputSTAR}"/RI_transcriptome"

if [ ! -d $outputSTAR ];then
	mkdir $outputSTAR
fi

if [ ! -d "$outputSTAR"/index ];then
	mkdir "$outputSTAR"/index
fi

if [ ! -d $nmd ];then
        mkdir $nmd
fi
if [ ! -d $ri ];then
        mkdir $ri
fi


#The following are the paths to the riboseq reads
#the two samples are merged
input_data="/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples/*.fastq"

#file with the 3'UTR background regiosn
three_primes="/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed"

#Create a Logfile for the alignments in the output directory
exec > >(tee -i $outputSTAR/AlignmentLogfile.txt)
exec 2>&1

#The following block calls the Bowtie_Align script, which creates a BOWTIE index for thetranscriptome and aligns
#Ribo-seq data against it before checking the overlap with the determined unique regions. 
#For further analysis a file with random regions from the 3' and 5' UTR
#is also created and used to determine background overlap

###everything in between should not be quoted, just to be faster#############################################
echo "Starting alignment against genome"


STAR --runThreadN 16 --runMode genomeGenerate --genomeDir "$outputSTAR"/index\
 --genomeFastaFiles /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa\
 --sjdbGTFfile /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf --sjdbOverhang 28




for i in $input_data; do
sample_name=$(basename $i _fastp.fastq)
source ./STAR_Align_genomic.sh  16 "$outputSTAR"/index $i\
 $nmd/${sample_name}_NMD $unique_region_dir/Unique_DNA_regions_genomic_NMD.bed $three_primes

source ./STAR_Align_genomic.sh 16 "$outputSTAR"/index $i\
 $ri/${sample_name}_RI $unique_region_dir/Unique_DNA_regions_genomic_RI.bed $three_primes

echo "===================       Sample $sample_name finished"

done





#Rscript -e "if (!requireNamespace('rmarkdown', quietly = TRUE)) install.packages('rmarkdown', repos='http://cran.us.r-project.org')"

#R -e "library(rmarkdown); rmarkdown::render(input = 'RiboSeqReportTranscriptomic_empirical_dist_dedup.Rmd', output_file = '$outputSTAR/Riboseq_report.pdf', params=list(args = c('$outputSTAR')))"

