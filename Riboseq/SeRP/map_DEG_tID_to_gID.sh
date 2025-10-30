#!/bin/bash

Bowtie2_out_dir=$1

eval "$(conda shell.bash hook)"
conda activate pygtftk

DEG_dir="${Bowtie2_out_dir}"/filtered/q10/DEGs

if [ ! -d $DEG_dir ]; then
        mkdir $DEG_dir
fi

enrichment_dir="/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/filtered/q10/enrichment_plots"


# obtain gene IDs of differentially expressed transcripts, CHX vs Mock
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_CHX_A2_mock_mRNA_0_5.txt  \
 $DEG_dir/DEGs_CHX_A2_mock_mRNA_0_5_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_CHX_B1_mock_mRNA_0_5.txt  \
 $DEG_dir/DEGs_CHX_B1_mock_mRNA_0_5_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_CHX_A2_mock_mRNA_0_5_downreg.txt  \
 $DEG_dir/DEGs_CHX_A2_mock_mRNA_0_5_downreg_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_CHX_B1_mock_mRNA_0_5_downreg.txt  \
 $DEG_dir/DEGs_CHX_B1_mock_mRNA_0_5_downreg_gene_IDs.txt 


# CHX vs Input
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_A2_In_CHX_mRNA_0_5.txt  \
 $DEG_dir/DEGs_A2_In_CHX_mRNA_0_5_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_B1_In_CHX_mRNA_0_5.txt  \
 $DEG_dir/DEGs_B1_In_CHX_mRNA_0_5_gene_IDs.txt 

 python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_A2_In_CHX_mRNA_0_5_downreg.txt  \
 $DEG_dir/DEGs_A2_In_CHX_mRNA_0_5_downreg_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_B1_In_CHX_mRNA_0_5_downreg.txt  \
 $DEG_dir/DEGs_B1_In_CHX_mRNA_0_5_downreg_gene_IDs.txt 


# Puro vs Input
python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_Puro_A2_mRNA_0_5.txt  \
 $DEG_dir/DEGs_Puro_A2_mRNA_0_5_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_Puro_B1_mRNA_0_5.txt  \
 $DEG_dir/DEGs_Puro_B1_mRNA_0_5_gene_IDs.txt 

 python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_Puro_A2_mRNA_0_5_downreg.txt  \
 $DEG_dir/DEGs_Puro_A2_mRNA_0_5_downreg_gene_IDs.txt 

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/DEG_analysis/map_tids_to_gids_gtf.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 $DEG_dir/DEGs_Puro_B1_mRNA_0_5_downreg.txt  \
 $DEG_dir/DEGs_Puro_B1_mRNA_0_5_downreg_gene_IDs.txt 



# Which genes are enriched in Mock and Input?
grep -Fxf $DEG_dir/DEGs_CHX_A2_mock_mRNA_0_5_gene_IDs.txt \
 $DEG_dir/DEGs_A2_In_CHX_mRNA_0_5_gene_IDs.txt \
 > $DEG_dir/DEGs_A2_CHX_0_5_enriched_over_Input_and_Mock_gene_IDs.txt


 grep -Fxf $DEG_dir/DEGs_CHX_B1_mock_mRNA_0_5_gene_IDs.txt \
 $DEG_dir/DEGs_B1_In_CHX_mRNA_0_5_gene_IDs.txt \
 > $DEG_dir/DEGs_B1_CHX__0_5enriched_over_Input_and_Mock_gene_IDs.txt



 # Enriched over both Input and Mock for 0.5
 grep -Fxf $DEG_dir/DEGs_CHX_A2_mock_mRNA_0_5.txt \
 $DEG_dir/DEGs_A2_In_CHX_mRNA_0_5.txt \
 > $DEG_dir/DEGs_A2_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs.txt


 grep -Fxf $DEG_dir/DEGs_CHX_B1_mock_mRNA_0_5.txt \
 $DEG_dir/DEGs_B1_In_CHX_mRNA_0_5.txt \
 > $DEG_dir/DEGs_B1_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs.txt

 grep -Fxf $DEG_dir/DEGs_A2_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs.txt \
 $DEG_dir/DEGs_B1_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs.txt \
> $DEG_dir/DEGs_both_A2_B1_CHX_0_5_Input_and_Mock_MANE_tIDs.txt


 # Enriched over both Input and Mock for 1
 grep -Fxf $DEG_dir/DEGs_CHX_A2_mock_mRNA_1.txt \
 $DEG_dir/DEGs_A2_In_CHX_mRNA_1.txt \
 > $DEG_dir/DEGs_A2_CHX_1_enriched_over_Input_and_Mock_MANE_tIDs.txt


 grep -Fxf $DEG_dir/DEGs_CHX_B1_mock_mRNA_1.txt \
 $DEG_dir/DEGs_B1_In_CHX_mRNA_1.txt \
 > $DEG_dir/DEGs_B1_CHX_1_enriched_over_Input_and_Mock_MANE_tIDs.txt

 grep -Fxf $DEG_dir/DEGs_A2_CHX_1_enriched_over_Input_and_Mock_MANE_tIDs.txt \
 $DEG_dir/DEGs_B1_CHX_1_enriched_over_Input_and_Mock_MANE_tIDs.txt \
> $DEG_dir/DEGs_both_A2_B1_CHX_1_Input_and_Mock_MANE_tIDs.txt


# Intersectoin with RIP seq targets
conda activate Riboseq
python map_gids_to_MANE_tids.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
 "${enrichment_dir}"/2025-06-04_RIP-Seq_hits.csv  \
 "${enrichment_dir}"/RIP_hits_MANE_tIDs.txt \
 "$DEG_dir"/RIP_hits_gids_to_MANE_tIDs.csv \
 multirow

 grep -Fxf $DEG_dir/DEGs_A2_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs.txt \
 $enrichment_dir/RIP_hits_MANE_tIDs.txt  \
 > $DEG_dir/DEGs_A2_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs_RIP_intersection.txt

grep -Fxf $DEG_dir/DEGs_B1_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs.txt \
 $enrichment_dir/RIP_hits_MANE_tIDs.txt  \
 > $DEG_dir/DEGs_B1_CHX_0_5_enriched_over_Input_and_Mock_MANE_tIDs_RIP_intersection.txt




DEG_files=($DEG_dir/*0_05.csv)
for csv in "${DEG_files[@]}"
do
python filter_DEG_csv_by_txt_index.py \
 ${csv} \
 $enrichment_dir/RIP_hits_MANE_tIDs.txt \
 $(dirname $csv)/$(basename $csv .csv)_RIP_histone_subset.csv
done