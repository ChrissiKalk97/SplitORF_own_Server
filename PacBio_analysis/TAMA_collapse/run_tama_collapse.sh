#!/bin/bash

#----- This script runs tama merge on 2 assemblies to merge them into one GTF file ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio


reference_gtf="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

bam_file_dir="/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment/HUVEC/pbmm2_align"

outdir_tama="/projects/splitorfs/work/PacBio/merged_bam_files/tama_collapse"

tama_tool_path="/home/ckalk/tools/tama"

bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/map_conditions/genome_mapping_cell_type.sh \
 -o "/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment" \
 -f "/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa" \
 -i "/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine" \
 -c "HUVEC"

bash tama_collapse.sh -b ${bam_file_dir} -c HUVEC -f ${genome_fasta}\
 -o ${outdir_tama} -r ${reference_gtf} -t ${tama_tool_path}