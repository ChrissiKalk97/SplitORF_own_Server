#!/bin/bash
#Help message:
usage="
Usage: ./Bowtie_Align_transcriptomic.sh [-options] numberOfThreads BowtieBaseName Reads.fastq out unique_regions.bed exonAnnotation genomicAnnotation

numberOfThreads 		Int setting the number of threads to use with BOWTIE2
BowtieBaseName 			BaseName to use for creating BOWTIE2 Index files, or if allready created, BaseName of the existing index files
Reads.fastq			Reads or transcripts that are to be aligned to the unique regions (can be fastq.gz)
out					Base name of the output files
unique_regions.bed		BED-file containing the unique_DNA_regions up for validation
transcripts.fa	fasta file with reference transcripts (NMD, RI etc.)

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

if [[ $# -ne 6 ]]; then #check for right number of arguments
  echo -e "${RED}
ERROR while executing the script!
Wrong number of arguments.${NC}"
  echo "Supplied arguments: $@" >&2 
  echo "$usage" >&2
  exit 1
fi


numberOfThreads=$1
StarIndex=$2
Riboreads=$3
bamfile=${4}.bam
unique_regions=$5
coordinates_3_prime=$6

echo "$Riboreads"

out_path=$(dirname $bamfile)


if [ ! -e  $out_path/genome_file.txt ]; then
    cut -f1,2 /projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa.fai >\
    $out_path/genome_chrom_ordering.txt
    chmod 777 $out_path/genome_chrom_ordering.txt
fi

if [ ! -e  $out_path/chromosomes_unique_regions.txt ]; then
    cut -f1 $unique_regions | sort | uniq > $out_path/chromosomes_unique_regions.txt
fi


sample=$(basename $bamfile .bam)
# # # align Riboreads against the genome
# STAR\
# --runThreadN $numberOfThreads\
# --alignEndsType EndToEnd\
# --outSAMstrandField intronMotif\
# --alignIntronMin 20\
# --alignIntronMax 1000000\
# --genomeDir $StarIndex\
# --readFilesIn $Riboreads\
# --twopassMode Basic\
# --outSAMattributes All\
# --outSAMtype BAM SortedByCoordinate\
# --outFileNamePrefix $out_path/$sample

# samtools view -@ $numberOfThreads -bo $bamfile $out_path/$(basename $bamfile .bam)Aligned.out.sam
bamfile=$out_path/${sample}Aligned.sortedByCoord.out.bam

filteredBamFile=$out_path/$(basename $bamfile Aligned.sortedByCoord.out.bam)_filtered.bam
# samtools view -F 256 -F 2048 -b $bamfile > $filteredBamFile
sortedBamFile=$out_path/$(basename $bamfile Aligned.sortedByCoord.out.bam)_sorted.bam
# samtools sort -o $sortedBamFile $filteredBamFile
# samtools index -@ 10 $sortedBamFile
bedfile=$out_path/$(basename $bamfile Aligned.sortedByCoord.out.bam).bed

present_chromosomes=$out_path/$(basename $bamfile Aligned.sortedByCoord.out.bam)_chromosomes.txt
# samtools view -H $sortedBamFile | grep '@SQ' | cut -f 2 | cut -d ':' -f 2  | sort | uniq > $present_chromosomes

# get the union of all chromosomes
# cat $present_chromosomes $out_path/chromosomes_unique_regions.txt | sort | uniq > temp_file
# mv temp_file $present_chromosomes


sorted_unique_regions=$(dirname $unique_regions)/$(basename $unique_regions .bed)_chrom_sorted.bed
if [ ! -e  $sorted_unique_regions ]; then
    sort -k1,1 -k2,2n $unique_regions > $sorted_unique_regions
fi



# echo "converting bam to bed"
# bedtools bamtobed -i $sortedBamFile -split > $bedfile


echo "intersecting with unique regions"
sortedBedfile=$out_path/$(basename $bedfile .bed)_intersect_counts_sorted.bed


# sort unique regions for intersect

#TEST#############

#############################################################
sorted_bedfile=$out_path/$(basename $bedfile .bed)_chrom_sort.bed
# sort -k1,1 -k2,2n $bedfile > $sorted_bedfile

# subset the genome bedfile to the present genomes
# grep -Fwf $present_chromosomes $out_path/genome_chrom_ordering.txt | sort -k1,1 -k2,2n > $out_path/genome_chrom_ordering_$(basename $bamfile Aligned.sortedByCoord.out.bam).txt
# -F: no regex, take chrs literally
# -w: match the whole word, e.g. 1 does not match 10, 11 etc but only 1
# -f: file input of the pattern that are searched for

# bedtools intersect\
#    -s\
#    -wao\
#    -a $sorted_unique_regions\
#    -b $sorted_bedfile\
#    -sorted\
#    -g $out_path/genome_chrom_ordering_$(basename $bamfile Aligned.sortedByCoord.out.bam).txt\
#    | sort -T /scratch/tmp/$USER -nr -k13,13\
#    > $sortedBedfile






for i in {1..20}; do
  randomfile=$out_path/Random_background_regions_${i}.bed
  sorted_randomfile="$out_path"/$(basename $randomfile .bed)_sorted_${i}.bed
  if [ ! -e  $randomfile ]; then
    python BackgroundRegions_bed_genomic_fix.py\
    $sorted_unique_regions\
    $coordinates_3_prime\
    $randomfile\
    $i

    # sort the randomfile with dummy entries
    echo $sorted_randomfile
    sort -T /scratch/tmp/"$USER" -k1,1 -k2,2n $randomfile > $sorted_randomfile
    rm $randomfile
  fi

  # echo $randomfile
  # echo $sorted_bedfile 
  # echo $sorted_randomfile

  # debug mode on
  # set -x 
  randomintersectfile=$out_path/$(basename $bamfile Aligned.sortedByCoord.out.bam)_${i}_random_intersect_counts.bed
  bedtools intersect\
    -s\
    -wao\
    -sorted\
    -g "$out_path/genome_chrom_ordering_$(basename $bamfile Aligned.sortedByCoord.out.bam).txt"\
    -a "$sorted_randomfile" \
    -b "$sorted_bedfile"\
    | sort -T /scratch/tmp/"$USER" -nr -k13,13\
    > $randomintersectfile
    # set +x 
    # debug mode off
done
# echo "Calculating random regions from 3 prime UTRs"




