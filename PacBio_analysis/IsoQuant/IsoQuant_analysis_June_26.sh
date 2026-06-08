#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate isoquant

genome_fasta_file="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
ensembl_gtf_filtered="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"

output_dir="/projects/splitorfs/work/PacBio/merged_bam_files/IsoQuant"

mkdir -p "$output_dir"


# assumpton: Preprocessing steps from 
# /home/ckalk/scripts/SplitORFs/PacBio_analysis/PacBio_analysis_29_05_26.sh
# are alaready run!

# might want to add this script as a subsript to the PacBio analysis



#################################################################################
# ------------------ Run IsoQuant                            ------------------ #
#################################################################################

isoquant -d pacbio_ccs \
 --polya_trimmed stranded \
 --sqanti_output \
 --yaml /home/ckalk/scripts/SplitORFs/PacBio_analysis/IsoQuant/long_reads.yaml  \
 --complete_genedb --genedb "${ensembl_gtf_filtered}" \
  --reference "$genome_fasta_file" --output "$output_dir"


 # --fl_data \
 # this option means that both ends of the reads are reliable, which is not the case