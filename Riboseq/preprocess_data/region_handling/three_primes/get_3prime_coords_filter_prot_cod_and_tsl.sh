#!/bin/bash -l
PATH="/home/fuchs/agschulz/kalk/miniforge3/bin:$PATH"
source /home/fuchs/agschulz/kalk/.bashrc
source /home/fuchs/agschulz/kalk/miniforge3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
source activate myenvname
python get_3prime_coords_filter_prot_cod_and_tsl.py\
 /scratch/fuchs/agschulz/kalk/star/filtered_Ensembl_reference/Ensembl_equality_and_TSL_filtered.gtf\
 /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/CDS_coords_110_no_filter.txt\
 /scratch/fuchs/agschulz/kalk/SplitORF/Riboseq/Input/three_primes_tsl1_refseq_prot_cod.bed

