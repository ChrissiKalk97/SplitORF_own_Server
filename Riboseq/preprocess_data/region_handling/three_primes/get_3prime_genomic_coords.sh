#!/bin/bash -l
eval "$(conda shell.bash hook)"
conda activate Riboseq
python get_3prime_genomic_coords.py\
 /projects/splitorfs/work/reference_files/Ensembl_equality_and_TSL_filtered.gtf\
 /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_prime_UTR_coords_genomic_Ensembl_110.txt\
 /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed

