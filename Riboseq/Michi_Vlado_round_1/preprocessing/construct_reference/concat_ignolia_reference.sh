#!/bin/bash

################################################################################
# Concatenate contaminants and mRNA according to Ignolia                       #
################################################################################

cat /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa \
  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta  \
  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa \
  /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta \
  > /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta



# python preprocessing/construct_reference/deduplicate_reference_fasta.py \
#  /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta \
#  /projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination_dedup.fasta
