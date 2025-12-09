#!/bin/bash
#Help message:
usage="
Usage: 

options...

where:
-h			show this help
"

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages


#available options for the programm
while getopts ':h' option; do
  case "$option" in
    h) echo "$usage"
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
  shift 2
done

if [[ $# -ne 6 ]]; then #check for right number of arguments
  echo -e "${RED}
  ERROR while executing the script!
  Wrong number of arguments.${NC}"
    echo "Supplied arguments: $@" >&2 
    echo "$usage" >&2
    exit 1
fi

################################################################################
# READ AND CHECK ARGUMENTS                                                     #
################################################################################

input_file=$1
unique_regions=$2
coordinates_3_prime=$3
outname=$4
random_region_path=$5
GENOME_FASTA=$6

echo "$input_file"
echo $unique_regions
echo $coordinates_3_prime
echo $outname
echo $random_region_path



sample=$(basename $outname)
out_path=$(dirname $outname)
echo $out_path


################################################################################
# GENERATE BED FILE IF NECESSARY                                              #
################################################################################

if [[ "$input_file" == *.bed ]]; then
  sorted_bedfile="$input_file"
else
    BED_FILE="$out_path/${sample}.bed"
    sorted_bedfile="$out_path/${sample}_chrom_sort.bed"
    bedtools bamtobed -i $input_file -split > $BED_FILE
    sort -k1,1 -k2,2n $BED_FILE > $sorted_bedfile
    rm $BED_FILE

    PRESENT_CHROMOSOMES="$out_path/$(basename $sorted_bedfile _chrom_sort.bed)_chromosomes.txt"
    samtools view -H $input_file | grep '@SQ' | cut -f 2 | cut -d ':' -f 2  | sort | uniq > $PRESENT_CHROMOSOMES
fi


################################################################################
# GENERATE REF FILES IF NECESSARY                                              #
################################################################################
if [ ! -s  $random_region_path/genome_file.txt ]; then
    cut -f1,2 ${GENOME_FASTA}.fai >\
    $random_region_path/genome_chrom_ordering.txt
    chmod 777 $random_region_path/genome_chrom_ordering.txt
fi


# chromosomes present in the unique regions
if [ ! -s  $random_region_path/chromosomes_unique_regions.txt ]; then
    cut -f1 $unique_regions | sort | uniq > $random_region_path/chromosomes_unique_regions.txt
fi

# sort unique regions according to chromosomes
sorted_unique_regions="$(dirname $unique_regions)/$(basename $unique_regions .bed)_chrom_sorted.bed"
if [ ! -s  $sorted_unique_regions ]; then
    sort -k1,1 -k2,2n $unique_regions > $sorted_unique_regions
fi


################################################################################
# INDEX, SORT, MAKE BED FROM GENOMIC ALIGNMENTS                                #
################################################################################
BED_DIR=$(dirname  $sorted_bedfile)

if [ ! -s  $BED_DIR/genome_chrom_ordering_$(basename $sorted_bedfile _chrom_sort.bed).txt ]; then
  PRESENT_CHROMOSOMES="$BED_DIR/$(basename $sorted_bedfile _chrom_sort.bed)_chromosomes.txt"
  grep -Fwf $PRESENT_CHROMOSOMES $random_region_path/genome_chrom_ordering.txt | sort -k1,1 -k2,2n > $BED_DIR/genome_chrom_ordering_$(basename $sorted_bedfile _chrom_sort.bed).txt
fi


echo "intersecting with unique regions"
intersectBedfile="${random_region_path}"/"${sample}"_intersect_counts_sorted.bed




################################################################################
# INTERSECT WITH UNIQUE REGIONS                                                #
################################################################################
bedtools intersect\
   -s\
   -wao\
   -a $sorted_unique_regions\
   -b $sorted_bedfile\
   -sorted\
   -g $BED_DIR/genome_chrom_ordering_$(basename $sorted_bedfile _chrom_sort.bed).txt\
   | sort -T /scratch/tmp/$USER -nr -k13,13\
   > $intersectBedfile




################################################################################
# INTERSECT WITH BACKGROUND REGIONS                                            #
################################################################################

for i in {1..20}; do
  randomfile="$random_region_path/Random_background_regions_${i}.bed"
  sorted_randomfile="$random_region_path/$(basename $randomfile .bed)_sorted_${i}.bed"
  if [ ! -s  $sorted_randomfile ]; then
    python /home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/BackgroundRegions_bed_genomic_fix.py\
    $sorted_unique_regions\
    $coordinates_3_prime\
    $randomfile\
    $i

    # sort the randomfile with dummy entries
    sort -T /scratch/tmp/"$USER" -k1,1 -k2,2n $randomfile > $sorted_randomfile
    echo $sorted_randomfile
    # rm $randomfile
  fi


  # debug mode on
  # set -x 
  randomintersectfile="${random_region_path}"/"${sample}"_${i}_random_intersect_counts.bed
  bedtools intersect\
    -s\
    -wao\
    -sorted\
    -g "$BED_DIR/genome_chrom_ordering_$(basename $sorted_bedfile _chrom_sort.bed).txt"\
    -a "$sorted_randomfile" \
    -b "$sorted_bedfile"\
    | sort -T /scratch/tmp/"$USER" -nr -k13,13\
    > $randomintersectfile


    # set +x 
    # debug mode off
done





