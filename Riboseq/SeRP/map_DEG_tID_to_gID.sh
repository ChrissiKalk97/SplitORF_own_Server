#!/bin/bash

Bowtie2_out_dir=$1

eval "$(conda shell.bash hook)"
conda activate pygtftk

if [ ! -d "${Bowtie2_out_dir}"/filtered/q10/DEGs ]; then
        mkdir "${Bowtie2_out_dir}"/filtered/q10/DEGs
fi


# obtain gene IDs of differentially expressed transcripts, CHX vs Mock
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_A2_mock_mRNA.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_A2_mock_mRNA_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_B1_mock_mRNA.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_B1_mock_mRNA_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_A2_mock_mRNA_downreg.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_A2_mock_mRNA_downreg_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_B1_mock_mRNA_downreg.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_B1_mock_mRNA_downreg_gene_IDs.txt 


# CHX vs Input
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_A2_In_CHX_mRNA.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_A2_In_CHX_mRNA_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_B1_In_CHX_mRNA.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_B1_In_CHX_mRNA_gene_IDs.txt 

 python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_A2_In_CHX_mRNA_downreg.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_A2_In_CHX_mRNA_downreg_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_B1_In_CHX_mRNA_downreg.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_B1_In_CHX_mRNA_downreg_gene_IDs.txt 


# Puro vs Input
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_Puro_A2_mRNA.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_Puro_A2_mRNA_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_Puro_B1_mRNA.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_Puro_B1_mRNA_gene_IDs.txt 

 python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_Puro_A2_mRNA_downreg.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_Puro_A2_mRNA_downreg_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_Puro_B1_mRNA_downreg.txt  \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_Puro_B1_mRNA_downreg_gene_IDs.txt 



# Which genes are enriched in Mock and Input?
grep -Fxf "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_A2_mock_mRNA_gene_IDs.txt \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_A2_In_CHX_mRNA_gene_IDs.txt \
 > "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_A2_CHX_enriched_over_Input_and_Mock_gene_IDs.txt


 grep -Fxf "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_B1_mock_mRNA_gene_IDs.txt \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_B1_In_CHX_mRNA_gene_IDs.txt \
 > "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_B1_CHX_enriched_over_Input_and_Mock_gene_IDs.txt



 # mRNA
 grep -Fxf "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_A2_mock_mRNA.txt \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_A2_In_CHX_mRNA.txt \
 > "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_A2_CHX_enriched_over_Input_and_Mock_MANE_tIDs.txt


 grep -Fxf "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_CHX_B1_mock_mRNA.txt \
 "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_B1_In_CHX_mRNA.txt \
 > "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_B1_CHX_enriched_over_Input_and_Mock_MANE_tIDs.txt


# Intersectoin with RIP seq targets
conda activate Riboseq
python map_gids_to_MANE_tids.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${Bowtie2_out_dir}"/filtered/q10/enrichment_plots/2025-06-04_RIP-Seq_hits.csv  \
 "${Bowtie2_out_dir}"/filtered/q10/enrichment_plots/RIP_hits_MANE_tIDs.txt \
 "${Bowtie2_out_dir}"/filtered/q10/enrichment_plots/RIP_hits_gids_to_MANE_tIDs.csv \
 multirow

 grep -Fxf "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_A2_CHX_enriched_over_Input_and_Mock_MANE_tIDs.txt \
 /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots/RIP_hits_MANE_tIDs.txt  \
 > /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_A2_CHX_enriched_over_Input_and_Mock_MANE_tIDs_RIP_intersection.txt

grep -Fxf "${Bowtie2_out_dir}"/filtered/q10/DEGs/DEGs_B1_CHX_enriched_over_Input_and_Mock_MANE_tIDs.txt \
 /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots/RIP_hits_MANE_tIDs.txt  \
 > /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/DEGs/DEGs_B1_CHX_enriched_over_Input_and_Mock_MANE_tIDs_RIP_intersection.txt



DEG_files=("${Bowtie2_out_dir}"/filtered/q10/DEGs/*0_05.csv)
for csv in "${DEG_files[@]}"
do
python filter_DEG_csv_by_txt_index.py \
 ${csv} \
 /projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots/RIP_hits_MANE_tIDs.txt \
 $(dirname $csv)/$(basename $csv .csv)_RIP_histone_subset.csv
done