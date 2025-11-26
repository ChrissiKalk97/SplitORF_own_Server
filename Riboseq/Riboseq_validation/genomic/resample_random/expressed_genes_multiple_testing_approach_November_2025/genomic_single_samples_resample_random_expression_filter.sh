#----- This script maps Ribo-seq data to the genome using the supplied annotation ----- #
# ----- then an intersection with unique regions in genome coords from the split-ORF pipeline ----- #
# ----- is performed as well as background regions of 3' UTRs, an empirical  ----- #
# ----- background distribution is used to determine which unique regions are  ----- #
# ----- and this is summarized in an Rmd report ----- #

#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq


output_star="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter"
unique_region_dir_nmd="/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-11.30.56_NMD_transcripts_correct_TSL_ref"
unique_region_dir_ri="/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-12.25.38_RI_transcripts_correct_TSL_ref"
ensembl_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
# The following are the paths to the riboseq reads
input_data=(/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples/*.fastq)
# file with the 3'UTR background regiosn
three_primes="/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed"

# please note: these are still aligned to Ens110
# would need to realign to respective index with TAMA GTF and also deduplicate
UMI_dedup_outdir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"
umi_dedup_upf10="/projects/splitorfs/work/UPF1_deletion/Output/alignment_genome/STAR/deduplicated"


if [ ! -d $output_star ];then
	mkdir $output_star
fi


if [ ! -d $output_star/NMD_genome ];then
	mkdir $output_star/NMD_genome
fi

if [ ! -d $output_star/RI_genome ];then
	mkdir $output_star/RI_genome
fi


# Create a Logfile for the alignments in the output directory
exec > >(tee -i $output_star/AlignmentLogfile.txt)
exec 2>&1

echo "STAR index genome"

if [ ! -d  "$output_star"/index ]; then
  source /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/STAR_Align_genomic_23_09_25.sh \
  -i 50 "$output_star"/index \
  $genome_fasta \
  $ensembl_gtf
fi

echo "Starting alignment against genome"


# STAR alignment for the samples that are "normal", no UMI deduplication etc
for i in "${input_data[@]}"; do
    sample_name=$(basename "$i" _fastp.fastq)
    
    if [ ! -e  "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_chrom_sort.bed ]; then
      echo $i
      mkdir "$output_star"/NMD_genome/${sample_name}

      /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/STAR_Align_genomic_23_09_25.sh \
      -a 16 "$output_star"/index $i \
      "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD \
      EndToEnd

      rm "$output_star"/NMD_genome/${sample_name}/${sample_name}*Aligned.sortedByCoord.out.bam
      rm "$output_star"/NMD_genome/${sample_name}/${sample_name}*_filtered.bam
      rm "$output_star"/NMD_genome/${sample_name}/*${sample_name}.bed

      echo "===================       Sample $sample_name mapped"
    fi

done



################################################################################
# NMD analysis                                                                 #
################################################################################
# Intersection for the "normal" Ribo-seq files
for folder in "$output_star"/NMD_genome/*/; do
    echo $folder
    bams=("$folder"/*NMD_sorted.bam)

    echo $bams

    if [[ -e "${bams[0]}" ]]; then
        i="${bams[0]}"        
        
        echo "here comes i"
        echo $i

        sample_name="$(basename "$i" _NMD_sorted.bam)"

        if [ ! -e  "$output_star"/NMD_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed ]; then

            htseq-count -f bam --secondary-alignments ignore \
              -c  "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv \
              --supplementary-alignments ignore $i $ensembl_gtf

            python filter_bed_file_for_expressed_genes.py \
              $three_primes \
              "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv

            python filter_bed_file_for_expressed_genes.py \
              $unique_region_dir_nmd/Unique_DNA_Regions_genomic.bed \
              "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv

            echo $sample_name
            ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh  \
                "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_chrom_sort.bed \
                "$output_star"/NMD_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed \
                "$output_star"/NMD_genome/${sample_name}/3_primes_genomic_${sample_name}.bed \
                "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD \
                "$output_star"/NMD_genome/${sample_name} \
                ${genome_fasta}

            echo "===================       Sample $sample_name intersected"
        fi
    fi

    

done

# the preprocessing and intersection for Vlado's protocol Ribo-seq data
# run this once in the first run, after that there are the chrom_sort.bed files
# already in the respective directory, and they are run in the second loop
for i in $UMI_dedup_outdir/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    if [ ! -e  "$output_star"/NMD_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed ]; then

        mkdir "$output_star"/NMD_genome/${sample_name}

        htseq-count -f bam --secondary-alignments ignore \
          -c  "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv \
          --supplementary-alignments ignore $i $ensembl_gtf

        python filter_bed_file_for_expressed_genes.py \
          $three_primes \
          "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv

        python filter_bed_file_for_expressed_genes.py \
          $unique_region_dir_nmd/Unique_DNA_Regions_genomic.bed \
          "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv

        echo $sample_name
        ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh  \
            $i \
            "$output_star"/NMD_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed \
            "$output_star"/NMD_genome/${sample_name}/3_primes_genomic_${sample_name}.bed \
            "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD \
            "$output_star"/NMD_genome/${sample_name} \
            ${genome_fasta}

        echo "===================       Sample $sample_name intersected"
    fi

done


# The preprocessing and intersection for the UPF1 deletion Ribo-seq data
for i in $umi_dedup_upf10/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    if [ ! -e  "$output_star"/NMD_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed ]; then
        mkdir "$output_star"/NMD_genome/${sample_name}

        htseq-count -f bam --secondary-alignments ignore \
          -c  "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv \
          --supplementary-alignments ignore $i $ensembl_gtf

        python filter_bed_file_for_expressed_genes.py \
          $three_primes \
          "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv

        python filter_bed_file_for_expressed_genes.py \
          $unique_region_dir_nmd/Unique_DNA_Regions_genomic.bed \
          "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD_htseq_counts.tsv

        echo $sample_name
        ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh  \
            $i \
            "$output_star"/NMD_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed \
            "$output_star"/NMD_genome/${sample_name}/3_primes_genomic_${sample_name}.bed \
            "$output_star"/NMD_genome/${sample_name}/${sample_name}_NMD \
            "$output_star"/NMD_genome/${sample_name} \
            ${genome_fasta}

        echo "===================       Sample $sample_name intersected"
    fi

done




export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
export MKL_ENABLE_INSTRUCTIONS=SSE4_2

Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos="http://cran.us.r-project.org")'

R -e 'library(rmarkdown); rmarkdown::render(input = "RiboSeqReportGenomic_iteration_update_expression_filter_multiple_test_correction.Rmd", output_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/Riboseq_report_NMD_expression_filter_multiple_test_correction.pdf", params=list(args = c("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome", "/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-11.30.56_NMD_transcripts_correct_TSL_ref", "NMD")))'





# ################################################################################
# # RI ANALYSIS                                                                  #
# ################################################################################
for folder in "$output_star"/NMD_genome/*/; do
    echo $folder
    bams=("$folder"/*NMD_sorted.bam)

    echo $bams

    if [[ -e "${bams[0]}" ]]; then
        i="${bams[0]}"        
        
        echo "here comes i"
        echo $i

        sample_name="$(basename "$i" _NMD_sorted.bam)"

        if [ ! -e  "$output_star"/RI_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed ]; then
            mkdir "$output_star"/RI_genome/${sample_name}

            htseq-count -f bam --secondary-alignments ignore \
              -c  "$output_star"/RI_genome/${sample_name}/${sample_name}_RI_htseq_counts.tsv \
              --supplementary-alignments ignore $i $ensembl_gtf

            python filter_bed_file_for_expressed_genes.py \
              $three_primes \
              "$output_star"/RI_genome/${sample_name}/${sample_name}_RI_htseq_counts.tsv

            python filter_bed_file_for_expressed_genes.py \
              $unique_region_dir_ri/Unique_DNA_Regions_genomic.bed \
              "$output_star"/RI_genome/${sample_name}/${sample_name}_RI_htseq_counts.tsv

            echo $sample_name
            ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh  \
                $i \
                "$output_star"/RI_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed \
                "$output_star"/RI_genome/${sample_name}/3_primes_genomic_${sample_name}.bed \
                "$output_star"/RI_genome/${sample_name}/${sample_name}_RI \
                "$output_star"/RI_genome/${sample_name} \
                ${genome_fasta}

            echo "===================       Sample $sample_name intersected"
        fi
    fi

    

done

# the preprocessing and intersection for Vlado's protocol Ribo-seq data
# run this once in the first run, after that there are the chrom_sort.bed files
# already in the respective directory, and they are run in the second loop
for i in $UMI_dedup_outdir/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    if [ ! -e  "$output_star"/RI_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed ]; then

        mkdir "$output_star"/RI_genome/${sample_name}

        htseq-count -f bam --secondary-alignments ignore \
          -c  "$output_star"/RI_genome/${sample_name}/${sample_name}_RI_htseq_counts.tsv \
          --supplementary-alignments ignore $i $ensembl_gtf

        python filter_bed_file_for_expressed_genes.py \
          $three_primes \
          "$output_star"/RI_genome/${sample_name}/${sample_name}_RI_htseq_counts.tsv

        python filter_bed_file_for_expressed_genes.py \
          $unique_region_dir_ri/Unique_DNA_Regions_genomic.bed \
          "$output_star"/RI_genome/${sample_name}/${sample_name}_RI_htseq_counts.tsv

        echo $sample_name
        ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh  \
            $i \
            "$output_star"/RI_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed \
            "$output_star"/RI_genome/${sample_name}/3_primes_genomic_${sample_name}.bed \
            "$output_star"/RI_genome/${sample_name}/${sample_name}_RI \
            "$output_star"/RI_genome/${sample_name} \
            ${genome_fasta}

        echo "===================       Sample $sample_name intersected"
    fi

done


# The preprocessing and intersection for the UPF1 deletion Ribo-seq data
for i in $umi_dedup_upf10/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    if [ ! -e  "$output_star"/RI_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed ]; then
        mkdir "$output_star"/RI_genome/${sample_name}

        htseq-count -f bam --secondary-alignments ignore \
          -c  "$output_star"/RI_genome/${sample_name}/${sample_name}_RI_htseq_counts.tsv \
          --supplementary-alignments ignore $i $ensembl_gtf

        python filter_bed_file_for_expressed_genes.py \
          $three_primes \
          "$output_star"/RI_genome/${sample_name}/${sample_name}_RI_htseq_counts.tsv

        python filter_bed_file_for_expressed_genes.py \
          $unique_region_dir_ri/Unique_DNA_Regions_genomic.bed \
          "$output_star"/RI_genome/${sample_name}/${sample_name}_RI_htseq_counts.tsv

        echo $sample_name
        ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh  \
            $i \
            "$output_star"/RI_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed \
            "$output_star"/RI_genome/${sample_name}/3_primes_genomic_${sample_name}.bed \
            "$output_star"/RI_genome/${sample_name}/${sample_name}_RI \
            "$output_star"/RI_genome/${sample_name} \
            ${genome_fasta}

        echo "===================       Sample $sample_name intersected"
    fi

done


export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
export MKL_ENABLE_INSTRUCTIONS=SSE4_2

Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos="http://cran.us.r-project.org")'

R -e 'library(rmarkdown); rmarkdown::render(input = "RiboSeqReportGenomic_iteration_update_expression_filter_multiple_test_correction.Rmd", output_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/Riboseq_report_RI_expression_filter_multiple_test_correction.pdf", params=list(args = c("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/RI_genome", "/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-11.30.56_NMD_transcripts_correct_TSL_ref", "RI")))'






# for i in "$output_star"/RI_genome/*_chrom_sort.bed; do
#     sample_name="$(basename "$i" _chrom_sort.bed)"

#     ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh \
#         $i \
#         $unique_region_dir_ri/Unique_DNA_Regions_genomic.bed \
#         $three_primes \
#         "$output_star"/RI_genome/${sample_name} \
#         "$output_star"/RI_genome \
#         ${genome_fasta}

#     echo "===================       Sample $sample_name intersected"

# done

# for i in "$UMI_dedup_outdir"/*_filtered.bam; do
#     sample_name=$(basename "$i" _dedup_filtered.bam)
    
#     echo $sample_name

#     if [[ "$sample_name" == uf_mueller* ]];then
#         sample_name="${sample_name:30}"
#     fi

#     echo $sample_name
#     ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh \
#         $i \
#         $unique_region_dir_ri/Unique_DNA_Regions_genomic.bed \
#         $three_primes \
#         "$output_star"/RI_genome/${sample_name}_RI \
#         "$output_star"/RI_genome \
#         ${genome_fasta}

#     echo "===================       Sample $sample_name intersected"

# done


# for i in "$umi_dedup_upf10"/*_filtered.bam; do
#     sample_name=$(basename "$i" _dedup_filtered.bam)
    
#     echo $sample_name

#     if [[ "$sample_name" == uf_mueller* ]];then
#         sample_name="${sample_name:30}"
#     fi

#     echo $sample_name
#     ./empirical_intersection_steps_expressed_genes_filter_18_11_25.sh \
#         $i \
#         $unique_region_dir_ri/Unique_DNA_Regions_genomic.bed \
#         $three_primes \
#         "$output_star"/RI_genome/${sample_name}_RI \
#         "$output_star"/RI_genome \
#         ${genome_fasta}

#     echo "===================       Sample $sample_name intersected"

# done



# export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
# export MKL_ENABLE_INSTRUCTIONS=SSE4_2

# Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos="http://cran.us.r-project.org")'

# R -e 'library(rmarkdown); rmarkdown::render(input = "RiboSeqReportGenomic_iteration_update_expression_filter_multiple_test_correction.Rmd", output_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/Riboseq_report_RI_expression_filter_multiple_test_correction.pdf", params=list(args = c("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/RI_genome", "/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-12.25.38_RI_transcripts_correct_TSL_ref", "RI")))'
