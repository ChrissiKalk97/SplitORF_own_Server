#!/bin/bash -l
eval "$(conda shell.bash hook)"
conda activate Riboseq
python get_3prime_genomic_coords.py\
 /projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf\
 /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_prime_UTR_coords_genomic_Ensembl_110.txt\
 /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed

sort -k1,1 -k2,2n /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed \
> /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_sorted.bed

awk 'BEGIN{OFS="\t"} {if($6=="1") $6="+"; else if($6=="-1") $6="-"; print}' \
 /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_sorted.bed \
  > /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_sorted_strand.bed

bedtools merge -i /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_sorted_strand.bed -s -c 4,5,6 -o collapse,min,distinct \
> /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_merged.bed

python number_genomic_regions.py \
 /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_merged.bed \
 /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_merged_numbered.bed