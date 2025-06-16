#!/bin/bash
#Help message:
usage="
Usage: 

options...

where:
-h			show this help
-i transcripts.fa	create new index files for the provided transcripts.fa"

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages


#available options for the programm
while getopts ':hi' option; do
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

if [[ $# -ne 5 ]]; then #check for right number of arguments
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

Ribobam=$1
unique_regions=$2
coordinates_3_prime=$3
outname=$4
random_region_path=$5

echo "$Ribobam"
echo $unique_regions
echo $coordinates_3_prime
echo $outname
echo $random_region_path



################################################################################
# GENERATE REF FILES IF NECESSARY                                              #
################################################################################

# chromosome ordering of reference genome, required for intersection
if [ ! -e  $random_region_path/genome_file.txt ]; then
    cut -f1,2 /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa.fai >\
    $random_region_path/genome_chrom_ordering.txt
    chmod 777 $random_region_path/genome_chrom_ordering.txt
fi

# chromosomes present in the unique regions
if [ ! -e  $random_region_path/chromosomes_unique_regions.txt ]; then
    cut -f1 $unique_regions | sort | uniq > $random_region_path/chromosomes_unique_regions.txt
fi

# sort unique regions according to chromosomes
sorted_unique_regions=$(dirname $unique_regions)/$(basename $unique_regions .bed)_chrom_sorted.bed
if [ ! -e  $sorted_unique_regions ]; then
    sort -k1,1 -k2,2n $unique_regions > $sorted_unique_regions
fi


################################################################################
# INDEX, SORT, MAKE BED FROM GENOMIC ALIGNMENTS                                #
################################################################################
sample=$(basename $Ribobam .bam)
out_path=$(dirname $outname)
echo $out_path

sortedBamFile=$out_path/$(basename $Ribobam .bam)_sorted.bam
# samtools sort -o $sortedBamFile $Ribobam
# samtools index -@ 10 $sortedBamFile

bedfile=$out_path/$(basename $Ribobam _dedup_filtered.bam).bed
if [ ! -e $bedfile ]; then
    echo "converting bam to bed"
    bedtools bamtobed -i $sortedBamFile -split > $bedfile
fi

present_chromosomes=$out_path/$(basename $Ribobam _dedup_filtered.bam)_chromosomes.txt
# samtools view -H $sortedBamFile | grep '@SQ' | cut -f 2 | cut -d ':' -f 2  | sort | uniq > $present_chromosomes



echo "intersecting with unique regions"
intersectBedfile="${random_region_path}"/"$(basename $outname)"_intersect_counts_sorted.bed


sorted_bedfile=$out_path/$(basename $Ribobam _dedup_filtered.bam)_chrom_sort.bed
if [ ! -e $out_path/$(basename $Ribobam _dedup_filtered.bam)_chrom_sort.bed ]; then
    sort -k1,1 -k2,2n $bedfile > $sorted_bedfile
fi


if [ ! -e $out_path/genome_chrom_ordering_$(basename $Ribobam _dedup_filtered.bam).txt ]; then
  # subset the genome bedfile to the present genomes
  grep -Fwf $present_chromosomes $random_region_path/genome_chrom_ordering.txt | sort -k1,1 -k2,2n > $out_path/genome_chrom_ordering_$(basename $Ribobam _dedup_filtered.bam).txt
  # -F: no regex, take chrs literally
  # -w: match the whole word, e.g. 1 does not match 10, 11 etc but only 1
  # -f: file input of the pattern that are searched for
fi



################################################################################
# INTERSECT WITH UNIQUE REGIONS                                                #
################################################################################
bedtools intersect\
   -s\
   -wao\
   -a $sorted_unique_regions\
   -b $sorted_bedfile\
   -sorted\
   -g $out_path/genome_chrom_ordering_$(basename $Ribobam _dedup_filtered.bam).txt\
   | sort -T /scratch/tmp/$USER -nr -k13,13\
   > $intersectBedfile




################################################################################
# INTERSECT WITH BACKGROUND REGIONS                                            #
################################################################################

for i in {1..20}; do
  randomfile=$random_region_path/Random_background_regions_${i}.bed
  sorted_randomfile=$random_region_path/$(basename $randomfile .bed)_sorted_${i}.bed
  if [ ! -e  $sorted_randomfile ]; then
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
  randomintersectfile="${random_region_path}"/"$(basename $outname)"_${i}_random_intersect_counts.bed
  bedtools intersect\
    -s\
    -wao\
    -sorted\
    -g "$out_path/genome_chrom_ordering_$(basename $Ribobam _dedup_filtered.bam).txt"\
    -a "$sorted_randomfile" \
    -b "$sorted_bedfile"\
    | sort -T /scratch/tmp/"$USER" -nr -k13,13\
    > $randomintersectfile


    # set +x 
    # debug mode off
done
# echo "Calculating random regions from 3 prime UTRs"




