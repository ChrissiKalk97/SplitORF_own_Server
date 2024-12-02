#!/bin/bash


#Help message:
usage="
Usage: ./STAR_Align_genomic.sh [-options] numberOfThreads STARBaseName Reads.fastq out unique_regions.bed three_primes.bed

numberOfThreads 		Int setting the number of threads to use with STAR
STARBaseName 			BaseName to use for creating STAR Index files, or if allready created, BaseName of the existing index files
Reads.fastq			Reads or transcripts that are to be aligned to the unique regions (can be fastq.gz), trimmed with fastp
out					Base name of the output files
unique_regions.bed		BED-file containing the unique_DNA_regions up for validation in genomic coordinates
three_primes.bed     BED-file containing the genomic background regions, 3' UTRs of protein coding tsl1 transcripts

where:
-h			show this help"
# -i transcripts.fa	create new index files for the provided transcripts.fa

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
  echo "$usage" >&2
  exit 1
fi


numberOfThreads=$1
StarIndex=$2
Riboreads=$3
bamfile=${4}.bam
unique_regions=$5
coordinates_3_prime=$6

out_path=$(dirname $bamfile)


if [ ! -e  $out_path/genome_file.txt ]; then
    cut -f1,2 /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa.fai >\
    $out_path/genome_chrom_ordering.txt
    chmod 777 $out_path/genome_chrom_ordering.txt
fi


# align Riboreads against the genome
STAR\
 --runThreadN $numberOfThreads\
 --alignEndsType EndToEnd\
 --outSAMstrandField intronMotif\
 --alignIntronMin 20\
 --alignIntronMax 1000000\
 --genomeDir $StarIndex\
 --readFilesIn $Riboreads\
 --twopassMode Basic\
 --outFileNamePrefix $out_path/$(basename $bamfile .bam)
#/scratch/fuchs/agschulz/kalk/star/reference_110_ribo
samtools view -@ $numberOfThreads -bo $bamfile $out_path/$(basename $bamfile .bam)Aligned.out.sam

filteredBamFile=$out_path/$(basename $bamfile .bam)_filtered.bam
samtools view -F 256 -F 2048 -b $bamfile > $filteredBamFile
sortedBamFile=$out_path/$(basename $bamfile .bam)_sorted.bam
samtools sort -o $sortedBamFile $filteredBamFile
samtools index -@ 10 $sortedBamFile
bedfile=$out_path/$(basename $bamfile .bam).bed

echo "converting bam to bed"
bedtools bamtobed -i $sortedBamFile -split > $bedfile


echo "intersecting with unique regions"
sortedBedfile=$out_path/$(basename $bedfile .bed)_intersect_counts_sorted.bed
sorted_unique_regions=$(dirname $unique_regions)/$(basename $unique_regions .bed)_chrom_sorted.bed

# sort unique regions for intersect
sort -k1,1 -k2,2n $unique_regions > $sorted_unique_regions
sorted_bedfile=$out_path/$(basename $bedfile .bed)_chrom_sort.bed
sort -k1,1 -k2,2n $bedfile > $sorted_bedfile


# intersect such that both entries are reported
# entry A is always reported with a null B feature, if no intersect
# file A the unique regions is a 6-bed, the Riboreads are a 12-bed
# the last col 19 gives the number of basepairs of overlap

bedtools intersect\
  -s\
  -wao\
  -F 0.33\
  -a $sorted_unique_regions\
  -b $sorted_bedfile\
  -sorted\
  -g $out_path/genome_chrom_ordering.txt\
  | sort -nr -k13,13\
  > $sortedBedfile



# intersectbedfilerelativesorted=$out_path/$(basename $bedfile .bed)_intersect_counts_relative_sorted.bed
# cat $sortedBedfile | awk -v OFS='\t' '{print $1,$2,$3,$4,$5, $5/($3-$2)}' | sort -nr -k 6 > $intersectbedfilerelativesorted 


echo "Calculating random regions from 3 prime UTRs"
randomfile=$out_path/$(basename $bamfile .bam)_random_background_regions.bed
python BackgroundRegions_bed_genomic.py\
 $sorted_unique_regions\
 $coordinates_3_prime\
 $randomfile

sorted_randomfile="$out_path"/$(basename $randomfile .bed)_sorted.bed
sort -k1,1 -k2,2n $randomfile > $sorted_randomfile

randomintersectfile=$out_path/$(basename $bamfile .bam)_random_intersect_counts.bed
bedtools intersect\
 -s\
  -wao\
  -F 0.33\
  -sorted\
  -g $out_path/genome_chrom_ordering.txt \
  -a $sorted_randomfile\
  -b $sorted_bedfile\
   | sort -nr -k13,13\
   > $randomintersectfile 
# randomintersectfilesorted=$out_path/$(basename $bamfile .bam)_random_intersect_counts_relative_sorted.bed
# cat $randomintersectfile | awk -v OFS='\t' '{print $1,$2,$3,$4, $4/($3-$2)}' | sort -n -r -k 5 > $randomintersectfilesorted
