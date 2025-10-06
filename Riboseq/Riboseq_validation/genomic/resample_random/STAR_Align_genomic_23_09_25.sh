#!/bin/bash
#Help message:
usage="
Usage: ./STAR_Align_genomic_23_09_25.sh [-options] [arguments]

where:
-h			show this help
-i NUMBER_OF_THREADS, STAR_INDEX, GENOME_FASTA and GENOME_GTF
-a NUMBER_OF_THREADS, STAR_INDEX, RIBO_READS, BAM_FILE, ALIGN_ENDS_TYPE"

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages


#available options for the programm
while getopts ':hi:a:' option; do
  case "$option" in
    h) 
        echo "$usage"
        exit 1
        ;;
    i)
        run_index=true
        if [ "$#" -lt 5 ]; then
            echo "Error: -i index requires 4 arguments: NUMBER_OF_THREADS, STAR_INDEX, GENOME_FASTA and GENOME_GTF" >&2
            exit 1
        fi

        # Assign them
        NUMBER_OF_THREADS=$2
        STAR_INDEX=$3
        GENOME_FASTA=$4
        GENOME_GTF=$5

        ;;

    a)
      run_align=true
      args=("$OPTARG" "${@:OPTIND:4}")   # get 4 args (1 from OPTARG, 3 more)
      if [ "${#args[@]}" -ne 5 ]; then
        echo "Error: -a requires 5 arguments: NUMBER_OF_THREADS, STAR_INDEX, RIBO_READS, BAM_FILE, ALIGN_ENDS_TYPE" >&2
        exit 1
      fi

      NUMBER_OF_THREADS=${args[0]}
      STAR_INDEX=${args[1]}
      RIBO_READS=${args[2]}
      BAM_NAME=${args[3]}
      ALIGN_ENDS_TYPE=${args[4]}
      ;;
   \?) 
        printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
  shift 2
done


if [ "$run_index" = true ]; then
        # index generation
        STAR --runThreadN 50 --runMode genomeGenerate --genomeDir "$STAR_INDEX" --genomeFastaFiles ${GENOME_FASTA}\
         --sjdbGTFfile ${GENOME_GTF} --sjdbOverhang 49
        # --sjdbOverhang 49: maxreadlength - 1: 50bp length limit for Ribo-seq data Vlado
fi

if [ "$run_align" = true ]; then

    echo "$RIBO_READS"

    OUT_PATH=$(dirname $BAM_NAME)


    if [ ! -e  $OUT_PATH/genome_file.txt ]; then
        cut -f1,2 ${GENOME_FASTA}.fai >\
        $OUT_PATH/genome_chrom_ordering.txt
        chmod 777 $OUT_PATH/genome_chrom_ordering.txt
    fi


    SAMPLE=$(basename $BAM_NAME)
    # align RIBO_READS against the genome
    STAR\
    --runThreadN $NUMBER_OF_THREADS\
    --alignEndsType ${ALIGN_ENDS_TYPE} \
    --outSAMstrandField intronMotif\
    --alignIntronMin 20\
    --alignIntronMax 1000000\
    --genomeDir $STAR_INDEX\
    --readFilesIn $RIBO_READS\
    --twopassMode Basic\
    --seedSearchStartLmax 20\
    --seedSearchStartLmaxOverLread 0.5\
    --outFilterMatchNminOverLread 0.9\
    --outSAMattributes All\
    --outSAMtype BAM SortedByCoordinate\
    --outFileNamePrefix ${BAM_NAME}

    BAM_FILE=${BAM_NAME}Aligned.sortedByCoord.out.bam

    FILTERED_BAM_FILE=${BAM_NAME}_filtered.bam
    samtools view -F 256 -F 2048 -q 10 -b $BAM_FILE > $FILTERED_BAM_FILE
    SORTED_BAM_FILE=${BAM_NAME}_sorted.bam
    samtools sort -o $SORTED_BAM_FILE $FILTERED_BAM_FILE
    samtools index -@ 10 $SORTED_BAM_FILE
    BED_FILE=${BAM_NAME}.bed

    PRESENT_CHROMOSOMES=${BAM_NAME}_chromosomes.txt
    samtools view -H $SORTED_BAM_FILE | grep '@SQ' | cut -f 2 | cut -d ':' -f 2  | sort | uniq > $PRESENT_CHROMOSOMES


    # echo "converting bam to bed"
    bedtools bamtobed -i $SORTED_BAM_FILE -split > $BED_FILE

    sorted_BED_FILE=$OUT_PATH/$(basename $BED_FILE .bed)_chrom_sort.bed

    sort -k1,1 -k2,2n $BED_FILE > $sorted_BED_FILE

    # subset the genome BED_FILE to the present genomes
    grep -Fwf $PRESENT_CHROMOSOMES $OUT_PATH/genome_chrom_ordering.txt | sort -k1,1 -k2,2n > $OUT_PATH/genome_chrom_ordering_${SAMPLE}.txt
    # -F: no regex, take chrs literally
    # -w: match the whole word, e.g. 1 does not match 10, 11 etc but only 1
    # -f: file input of the pattern that are searched for
fi









