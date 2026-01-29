#----- This script maps Ribo-seq data to the genome using the supplied annotation ----- #
# ----- then an intersection with unique regions in genome coords from the split-ORF pipeline ----- #
# ----- is performed as well as background regions of 3' UTRs, an empirical  ----- #
# ----- background distribution is used to determine which unique regions are  ----- #
# ----- and this is summarized in an Rmd report ----- #

#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq


output_star="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter"
unique_region_dir_ri="/home/ckalk/tools/SplitORF_pipeline/Output/run_13.01.2026-12.48.25_RI_CDS_subtraction_minAlignLength_15"
ensembl_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
# The following are the paths to the riboseq reads
input_data=(/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples/*.fastq)
# file with the 3'UTR background regiosn
three_primes="/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_merged_numbered.bed"

# please note: these are still aligned to Ens110
# would need to realign to respective index with TAMA GTF and also deduplicate
UMI_dedup_outdir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"
umi_dedup_upf10="/projects/splitorfs/work/UPF1_deletion/Output/alignment_genome/STAR/deduplicated"

cds_coordinates="/projects/splitorfs/work/Riboseq/data/region_input/genomic/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered.bed"


if [ ! -d $output_star ];then
	mkdir $output_star
fi

input_name="NMD"
file_dir="${output_star}/"${input_name}"_genome"

region_type="RI"




# keep the NMD genome: will place alignments here
if [ ! -d "${file_dir}" ];then
	mkdir "${file_dir}"
fi



# getting the 3' UTR regions
bash /home/ckalk/scripts/SplitORFs/Riboseq/preprocess_data/region_handling/three_primes/get_3prime_genomic_coords.sh



# Create a Logfile for the alignments in the output directory
exec > >(tee -i $output_star/AlignmentLogfile.txt)
exec 2>&1

echo "STAR index genome"

if [ ! -d  "$output_star"/index ]; then
  source /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/STAR_Align_genomic_23_09_25.sh \
  -i 50 "$output_star"/index \
  "$genome_fasta" \
  "$ensembl_gtf"
fi

echo "Starting alignment against genome"


# STAR alignment for the samples that are "normal", no UMI deduplication etc
for i in "${input_data[@]}"; do
    sample_name=$(basename "$i" _fastp.fastq)
    
    if [[ ! -e  ""${file_dir}"/${sample_name}/${sample_name}_"${input_name}"_chrom_sort.bed" ]]; then
      echo $i
      mkdir "${file_dir}"/${sample_name}

      /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/STAR_Align_genomic_23_09_25.sh \
      -a 16 "$output_star"/index $i \
      "${file_dir}"/${sample_name}/${sample_name}_"${input_name}" \
      EndToEnd \
      "$genome_fasta"

      rm "${file_dir}"/${sample_name}/${sample_name}*Aligned.sortedByCoord.out.bam
      rm "${file_dir}"/${sample_name}/${sample_name}*_filtered.bam

      echo "===================       Sample $sample_name mapped"
    fi

done

# Preparing files for Picard
if [[ ! -e $ensembl_gtf.refflat ]]; then
    gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons $ensembl_gtf /dev/stdout | \
        awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > $ensembl_gtf.refflat
fi



################################################################################
# NMD analysis                                                                 #
################################################################################
# Intersection for the "normal" Ribo-seq files
for folder in "${file_dir}"/*/; do
    echo $folder
    bams=("$folder"/*"${input_name}"_sorted.bam)

    echo $bams

    if [[ -e "${bams[0]}" ]]; then
        i="${bams[0]}"   


        sample_name="$(basename "$i" _"${input_name}"_sorted.bam)"

        echo $sample_name

        bash filter_intersection_pipeline_region_type.sh \
        -b $i \
        -c $cds_coordinates \
        -e $ensembl_gtf\
        -g $genome_fasta \
        -i "${file_dir}/${sample_name}/${sample_name}_${input_name}_chrom_sort.bed" \
        -n "${file_dir}/${sample_name}/${sample_name}_${input_name}" \
        -o $output_star \
        -s $sample_name \
        -r "${region_type}" \
        -t $three_primes \
        -u $unique_region_dir_ri \
        -d
    fi

    if [[ ! -e "${output_star}/"${input_name}"_genome/RiboseQC/ORFquant/RiboseQC/${sample_name}_04_08_25.html" ]]; then
        bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/Riboseq_pipeline.sh \
        -i "${output_star}/"${input_name}"_genome/${sample_name}" \
        -o "${file_dir}/RiboseQC" \
        -g "${ensembl_gtf}" \
        -q
    fi
done






# the preprocessing and intersection for Vlado's protocol Ribo-seq data
# run this once in the first run, after that there are the chrom_sort.bed files
# already in the respective directory, and they are run in the second loop
for i in "${UMI_dedup_outdir}"/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    if [ ! -e  "$output_star"/"${region_type}"_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed ]; then

        if [ ! -d "$output_star"/"${region_type}"_genome/${sample_name} ];then
            mkdir "$output_star"/"${region_type}"_genome/${sample_name}
        fi


        bash filter_intersection_pipeline_region_type.sh \
        -b $i \
        -c $cds_coordinates \
        -e $ensembl_gtf\
        -g $genome_fasta \
        -i $i \
        -n "${file_dir}/${sample_name}/${sample_name}_${input_name}" \
        -o $output_star \
        -s $sample_name \
        -r "${region_type}" \
        -t $three_primes \
        -u $unique_region_dir_ri

    fi

done

if [[ ! -e "${file_dir}/RiboseQC/ORFquant/RiboseQC/${sample_name}_dedup_filtered.bam_04_08_25.html" ]]; then
    bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/Riboseq_pipeline.sh \
        -i  "${UMI_dedup_outdir}" \
        -o "${file_dir}/RiboseQC" \
        -g "${ensembl_gtf}" \
        -q
fi


# The preprocessing and intersection for the UPF1 deletion Ribo-seq data
for i in "${umi_dedup_upf10}"/*_filtered.bam; do
    sample_name=$(basename "$i" _dedup_filtered.bam)
    
    echo $sample_name

    if [[ "$sample_name" == uf_mueller* ]];then
        sample_name="${sample_name:30}"
    fi

    if [ ! -e  "$output_star"/"${region_type}"_genome/${sample_name}/Unique_DNA_Regions_genomic_${sample_name}.bed ]; then
        if [ ! -d "$output_star"/"${region_type}"_genome/${sample_name} ];then
            mkdir "$output_star"/"${region_type}"_genome/${sample_name}
        fi

        bash filter_intersection_pipeline_region_type.sh \
        -b $i \
        -c $cds_coordinates \
        -e $ensembl_gtf \
        -g $genome_fasta \
        -i $i \
        -n "${file_dir}/${sample_name}/${sample_name}_${input_name}" \
        -o $output_star \
        -s $sample_name \
        -r "${region_type}" \
        -t $three_primes \
        -u $unique_region_dir_ri
        
    fi

done

if [[ ! -e "${file_dir}/RiboseQC/ORFquant/RiboseQC/${sample_name}_04_08_25.html" ]]; then
    bash /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/Riboseq_pipeline.sh \
        -i  "${umi_dedup_upf10}" \
        -o "${file_dir}/RiboseQC" \
        -g "${ensembl_gtf}" \
        -q
fi


if [ ! -e "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/Riboseq_report_NMD_3prime_CDS_quantile_filter_20bp_windows_13_01_26.pdf" ]; then
    export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
    export MKL_ENABLE_INSTRUCTIONS=SSE4_2

    Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos="http://cran.us.r-project.org")'

    R -e 'library(rmarkdown); rmarkdown::render(input = "RiboSeqReportGenomic_iteration_update_expression_filter_multiple_test_correction.Rmd", output_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/Riboseq_report_RI_filter_sample_7mio_threshold_19_01_26.pdf", params=list(args = c("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/RI_genome", "/home/ckalk/tools/SplitORF_pipeline/Output/run_13.01.2026-12.48.25_RI_CDS_subtraction_minAlignLength_15", "RI")))'
fi

# get a summary of the transcripts that have 2 URs or more validated
python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/expressed_genes_multiple_testing_approach_November_2025/downstream_analysis_validated_URs/summarize_transcripts_with_2_val_regions_across_samples.py \
    --out_name '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/RI_genome/transcripts_with_2_regions_summarized.csv' \
    --ribo_coverage_path '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/RI_genome' \
    --samples_of_interest_string 'SRR10,SRR85,HCT'

