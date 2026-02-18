#!/bin/bash -l

# This script intersects each genome aligned bam file with the unique
# regions in genomic coordinats
# For the background estiamtion a file with random regions from the 3' and 5' UTR
# is also created and used to determine background overlap
# An RMD report of the significance of the results is created in the end


eval "$(conda shell.bash hook)"
source activate Riboseq


outputSTAR="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"
nmd=${outputSTAR}"/empirical_Riboseq_validation/NMD_genome"
unique_region_dir="/projects/splitorfs/work/Riboseq/data/region_input/genomic"
ri=${outputSTAR}"/empirical_Riboseq_validation/RI_genome"
random_region_path_NMD="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/NMD_genome"
random_region_path_RI="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/RI_genome"


if [ ! -d $outputSTAR ];then
	mkdir $outputSTAR
fi

if [ ! -d "${outputSTAR}"/empirical_Riboseq_validation ];then
	mkdir "${outputSTAR}"/empirical_Riboseq_validation
fi

if [ ! -d $nmd ];then
        mkdir $nmd
fi
if [ ! -d $ri ];then
        mkdir $ri
fi


 #The following are the paths to the riboseq reads
#the two samples are merged
input_data="${outputSTAR}"/*filtered.bam

#file with the 3'UTR background regiosn
three_primes="/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed"




###everything in between should not be quoted, just to be faster#############################################
echo "Starting empirical Riboseq validation"



for Ribobam in $input_data; do
        echo $Ribobam
        sample_name=$(basename "$Ribobam" _dedup_filtered.bam)

        bash ./Riboseq_SO_empirical_intersection/empirical_intersection_random_resample.sh \
        "$Ribobam" \
        "$unique_region_dir"/Unique_DNA_regions_genomic_NMD_16_12_24.bed \
        "$three_primes" \
        "$nmd"/"${sample_name}"_NMD \
        "$random_region_path_NMD"

        echo "===================       Sample $sample_name NMD finished"

done


for Ribobam in $input_data; do
        echo $Ribobam
        sample_name=$(basename "$Ribobam" _dedup_filtered.bam)

        bash ./Riboseq_SO_empirical_intersection/empirical_intersection_random_resample.sh \
        "$Ribobam" \
        "$unique_region_dir"/Unique_DNA_regions_genomic_RI_16_12_24.bed \
        "$three_primes" \
        "$ri"/"${sample_name}"_RI \
        "$random_region_path_RI"

        echo "===================       Sample $sample_name RI finished"

done





export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
export MKL_ENABLE_INSTRUCTIONS=SSE4_2

Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos="http://cran.us.r-project.org")'

cd /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random
R -e 'library(rmarkdown); rmarkdown::render(input = "RiboSeqReportGenomic_iteration.Rmd", output_file = "/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated/empirical_Riboseq_validation/Riboseq_emp_val_report_hypo_q10_16_06_15.pdf", params=list(args = c("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/", "/projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results")))'