#----- This script takes the genomic unique regions from the Split-ORf pipeline ----- #
#----- and checks for their overlap with CDS regions (protein coding, TSL and RefSeq fitlered) ----- #
#----- It also removes the overlapping regions (bedtools subtract) to enable a "clean version" ----- #


#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate pygtftk

cds_coordinate_bed="/projects/splitorfs/work/Riboseq/data/region_input/genomic/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered.bed"

if [ ! -e ${cds_coordinate_bed} ]; then
    python Ens_CDS_coord_bed_gtf_filtered.py \
        /projects/splitorfs/work/Riboseq/data/region_input/genomic/Ens_110_CDS_coordinates_genomic_all.txt \
        /projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf \
        ${cds_coordinate_bed}
fi

conda activate Riboseq


unique_region_dir_nmd="/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-11.30.56_NMD_transcripts_correct_TSL_ref"
unique_region_dir_ri="/home/ckalk/tools/SplitORF_pipeline/Output/run_30.09.2025-12.25.38_RI_transcripts_correct_TSL_ref"

bedtools intersect -s -wo -a ${unique_region_dir_nmd}/Unique_DNA_Regions_genomic.bed -b ${cds_coordinate_bed} > ${unique_region_dir_nmd}/Unique_DNA_Regions_genomic_CDS_intersection.bed

bedtools subtract -s -a ${unique_region_dir_nmd}/Unique_DNA_Regions_genomic.bed -b ${cds_coordinate_bed} > ${unique_region_dir_nmd}/Unique_DNA_Regions_genomic_CDS_subtraction.bed

python count_CDS_overlapping_unique_regions.py \
 ${unique_region_dir_nmd}/Unique_DNA_Regions_genomic_CDS_intersection.bed


 bedtools intersect -s -wo -a ${unique_region_dir_ri}/Unique_DNA_Regions_genomic.bed -b ${cds_coordinate_bed} > ${unique_region_dir_ri}/Unique_DNA_Regions_genomic_CDS_intersection.bed

bedtools subtract -s -a ${unique_region_dir_ri}/Unique_DNA_Regions_genomic.bed -b ${cds_coordinate_bed} > ${unique_region_dir_ri}/Unique_DNA_Regions_genomic_CDS_subtraction.bed

python count_CDS_overlapping_unique_regions.py \
 ${unique_region_dir_ri}/Unique_DNA_Regions_genomic_CDS_intersection.bed

