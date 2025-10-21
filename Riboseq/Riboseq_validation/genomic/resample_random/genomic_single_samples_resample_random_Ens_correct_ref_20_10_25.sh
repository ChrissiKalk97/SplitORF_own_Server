#----- This script maps Ribo-seq data to the genome using the supplied annotation ----- #
# ----- then an intersection with unique regions in genome coords from the split-ORF pipeline ----- #
# ----- is performed as well as background regions of 3' UTRs, an empirical  ----- #
# ----- background distribution is used to determine which unique regions are  ----- #
# ----- and this is summarized in an Rmd report ----- #

#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq


output_star="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10"
unique_region_dir_nmd="/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-11.30.56_NMD_transcripts_correct_TSL_ref"
unique_region_dir_ri="/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-12.25.38_RI_transcripts_correct_TSL_ref"
ensembl_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
# The following are the paths to the riboseq reads
# INPUT_DATA=(/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples/*.fastq)
# file with the 3'UTR background regiosn
three_primes="/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed"

# please note: these are still aligned to Ens110
# would need to realign to respective index with TAMA GTF and also deduplicate
UMI_dedup_outdir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"
umi_dedup_upf10="/projects/splitorfs/work/UPF1_deletion/Output/alignment_genome/STAR/deduplicated"


if [ ! -d $output_star ];then
	mkdir $output_star
fi



# Create a Logfile for the alignments in the output directory
exec > >(tee -i $output_star/AlignmentLogfile.txt)
exec 2>&1

# echo "STAR index genome"

# source ./STAR_Align_genomic_23_09_25.sh -i 50 "$output_star"/index \
#  $genome_fasta \
#  $ensembl_gtf

# echo "Starting alignment against genome"

# for i in "${INPUT_DATA[@]}"; do
#     sample_name=$(basename "$i" _fastp.fastq)
#     echo $i

#     ./STAR_Align_genomic_23_09_25.sh -a 16 "$output_star"/index $i \
#     "$output_star"/NMD_genome/${sample_name} \
#      EndToEnd

#     rm "$output_star"/NMD_genome/${sample_name}*Aligned.sortedByCoord.out.bam
#     rm "$output_star"/NMD_genome/${sample_name}*_filtered.bam
#     rm "$output_star"/NMD_genome/*${sample_name}.bed

#     echo "===================       Sample $sample_name mapped"

# done



################################################################################
# NMD analysis                                                                 #
################################################################################

for i in $output_star/NMD_genome/*_chrom_sort.bed; do
    sample_name="$(basename "$i" _chrom_sort.bed)"

    ./empirical_intersection_steps_23_09_25.sh \
        $i \
        $unique_region_dir_nmd/Unique_DNA_Regions_genomic.bed \
        $three_primes \
        "$output_star"/NMD_genome/${sample_name} \
        "$output_star"/NMD_genome \
        ${genome_fasta}

    echo "===================       Sample $sample_name intersected"

done


# run this once in the first run, after that there are the chrom_sort.bed files
# already in the respective directory, and they are run in the second loop
for i in $UMI_dedup_outdir/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    echo $sample_name
    ./empirical_intersection_steps_23_09_25.sh \
        $i \
        $unique_region_dir_nmd/Unique_DNA_Regions_genomic.bed \
        $three_primes \
        "$output_star"/NMD_genome/${sample_name}_NMD \
        "$output_star"/NMD_genome \
        ${genome_fasta}

    echo "===================       Sample $sample_name intersected"

done


for i in $umi_dedup_upf10/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    echo $sample_name
    ./empirical_intersection_steps_23_09_25.sh \
        $i \
        $unique_region_dir_nmd/Unique_DNA_Regions_genomic.bed \
        $three_primes \
        "$output_star"/NMD_genome/${sample_name}_NMD \
        "$output_star"/NMD_genome \
        ${genome_fasta}

    echo "===================       Sample $sample_name intersected"

done




export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
export MKL_ENABLE_INSTRUCTIONS=SSE4_2

Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos="http://cran.us.r-project.org")'

R -e 'library(rmarkdown); rmarkdown::render(input = "RiboSeqReportGenomic_iteration_update_single_23_09_25.Rmd", output_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/Riboseq_report_NMD_Ens_correct_ref_20_10_25.pdf", params=list(args = c("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/NMD_genome", "/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-11.30.56_NMD_transcripts_correct_TSL_ref", "NMD")))'





################################################################################
# RI ANALYSIS                                                                  #
################################################################################

for i in "$output_star"/RI_genome/*_chrom_sort.bed; do
    sample_name="$(basename "$i" _chrom_sort.bed)"

    ./empirical_intersection_steps_23_09_25.sh \
        $i \
        $unique_region_dir_ri/Unique_DNA_Regions_genomic.bed \
        $three_primes \
        "$output_star"/RI_genome/${sample_name} \
        "$output_star"/RI_genome \
        ${genome_fasta}

    echo "===================       Sample $sample_name intersected"

done

for i in "$UMI_dedup_outdir"/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    echo $sample_name
    ./empirical_intersection_steps_23_09_25.sh \
        $i \
        $unique_region_dir_ri/Unique_DNA_Regions_genomic.bed \
        $three_primes \
        "$output_star"/RI_genome/${sample_name}_RI \
        "$output_star"/RI_genome \
        ${genome_fasta}

    echo "===================       Sample $sample_name intersected"

done


for i in "$umi_dedup_upf10"/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    echo $sample_name
    ./empirical_intersection_steps_23_09_25.sh \
        $i \
        $unique_region_dir_ri/Unique_DNA_Regions_genomic.bed \
        $three_primes \
        "$output_star"/RI_genome/${sample_name}_RI \
        "$output_star"/RI_genome \
        ${genome_fasta}

    echo "===================       Sample $sample_name intersected"

done



export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
export MKL_ENABLE_INSTRUCTIONS=SSE4_2

Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos="http://cran.us.r-project.org")'

R -e 'library(rmarkdown); rmarkdown::render(input = "RiboSeqReportGenomic_iteration_update_single_23_09_25.Rmd", output_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/Riboseq_report_RI_Ens_correct_ref_20_10_25.pdf", params=list(args = c("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/RI_genome", "/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-12.25.38_RI_transcripts_correct_TSL_ref", "RI")))'
