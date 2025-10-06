#----- This script maps Ribo-seq data to the genome using the supplied annotation ----- #
# ----- then an intersection with unique regions in genome coords from the split-ORF pipeline ----- #
# ----- is performed as well as background regions of 3' UTRs, an empirical  ----- #
# ----- background distribution is used to determine which unique regions are  ----- #
# ----- and this is summarized in an Rmd report ----- #

#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq


OUTPUT_STAR="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_huvec_tama_19_09_25"
UNIQUE_REGION_DIR="/projects/splitorfs/work/Masspec/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25"
TAMA_GTF="/projects/splitorfs/work/PacBio/merged_bam_files/compare_mando_stringtie/tama/HUVEC/HUVEC_merged_tama_gene_id.gtf"
GENOME_FASTA="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
# The following are the paths to the riboseq reads
INPUT_DATA=(/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples/*.fastq)
# file with the 3'UTR background regiosn
THREE_PRIMES="/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed"

# please note: these are still aligned to Ens110
# would need to realign to respective index with TAMA GTF and also deduplicate
UMI_dedup_outdir="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome_huvec_tama/STAR/only_R1/deduplicated"


if [ ! -d $OUTPUT_STAR ];then
	mkdir $OUTPUT_STAR
fi



# Create a Logfile for the alignments in the output directory
# exec > >(tee -i $OUTPUT_STAR/AlignmentLogfile.txt)
# exec 2>&1

# echo "STAR index genome"

# source ./STAR_Align_genomic_23_09_25.sh -i 50 "$OUTPUT_STAR"/index \
#  $GENOME_FASTA \
#  $TAMA_GTF

# echo "Starting alignment against genome"

# for i in "${INPUT_DATA[@]}"; do
#     sample_name=$(basename "$i" _fastp.fastq)
#     echo $i

#     ./STAR_Align_genomic_23_09_25.sh -a 16 "$OUTPUT_STAR"/index $i \
#     "$OUTPUT_STAR"/${sample_name} \
#      EndToEnd

#     rm "$OUTPUT_STAR"/${sample_name}*Aligned.sortedByCoord.out.bam
#     rm "$OUTPUT_STAR"/${sample_name}*_filtered.bam
#     rm "$OUTPUT_STAR"/*${sample_name}.bed

#     echo "===================       Sample $sample_name mapped"

# done

# # run this once in the first run, after that there are the chrom_sort.bed files
# # already in the respective directory, and they are run in the second loop
# for i in $UMI_dedup_outdir/*_filtered.bam; do
#     sample_name=$(basename "$i" _dedup_filtered.bam)
    
#     echo $sample_name

#     if [[ "$sample_name" == uf_mueller* ]];then
#         sample_name="${sample_name:30}"
#     fi

#     echo $sample_name
#    ./empirical_intersection_steps_23_09_25.sh \
#         $i \
#         $UNIQUE_REGION_DIR/huvec_tama_unique_peptides_masspec_genomic.bed \
#         $THREE_PRIMES \
#         "$OUTPUT_STAR"/${sample_name}_HUVEC_TAMA \
#         "$OUTPUT_STAR" \
#         ${GENOME_FASTA}

#     echo "===================       Sample $sample_name intersected"

# done

# ################################################################################
# # PREPARE BED FILES
# ################################################################################
# if [ ! -s  $UNIQUE_REGION_DIR/huvec_tama_unique_peptides_masspec_genomic.bed ]; then
# 	cat $UNIQUE_REGION_DIR/*genomic_coordinates.bed | sort -k1,1 -k2,2n > $UNIQUE_REGION_DIR/huvec_tama_unique_peptides_masspec_genomic.bed
# fi

# for i in $OUTPUT_STAR/*_chrom_sort.bed; do
#     sample_name="$(basename "$i" _chrom_sort.bed)"

#     ./empirical_intersection_steps_23_09_25.sh \
#         $i \
#         $UNIQUE_REGION_DIR/huvec_tama_unique_peptides_masspec_genomic.bed \
#         $THREE_PRIMES \
#         "$OUTPUT_STAR"/HUVEC_peptides/${sample_name}_HUVEC_PEPTIDES \
#         "$OUTPUT_STAR"/HUVEC_peptides \
#         ${GENOME_FASTA}

#     echo "===================       Sample $sample_name intersected"

# done



# export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64:$LD_LIBRARY_PATH
# export MKL_ENABLE_INSTRUCTIONS=SSE4_2

# Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos="http://cran.us.r-project.org")'

# R -e 'library(rmarkdown); rmarkdown::render(input = "RiboSeqReportGenomic_iteration_update_single_23_09_25.Rmd", output_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_huvec_tama_19_09_25/HUVEC_peptides/Riboseq_report_TAMA_HUVEC_peptides.pdf", params=list(args = c("/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_huvec_tama_19_09_25/HUVEC_peptides", "/home/ckalk/tools/SplitORF_pipeline/Output/run_12.09.2025-17.51.04_HUVEC_tama_merged", "HUVEC_PEPTIDES")))'

# Rscript plotting/Upsetplot_SO_upset_by_type_ISMB_talk_01_07_25.R $nmd $ri TRUE
# Rscript plotting/Upsetplot_SO_upset_by_type_ISMB_talk_01_07_25.R $nmd $ri FALSE


################################################################################
# GET GENE IDS FOR GO AND STRINGS OF PROTEINS WITH RIBOSEQ COVERAGE
################################################################################

python unique_peptide_helper_scripts/get_validated_protein_genes_masspec.py \
    /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_huvec_tama_19_09_25/HUVEC_peptides \
    /projects/splitorfs/work/Masspec/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/huvec_validated_SO_protein_original_Ids_with_assembly.csv \
    /projects/splitorfs/work/Masspec/New_MS_run_19_09_25_tama_assembly_SOs/analysis_results_with_ref_19_09_25/huvec_tama_unique_peptides_masspec_genomic_chrom_sorted.bed
