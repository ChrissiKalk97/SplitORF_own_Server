#!/bin/bash


# Help message:
usage="
Usage: ./genome_alignment_star.sh [-options] [arguments]

where:
-h			show this help
"

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages


#available options for the programm
while getopts 'a:e:f:g:him:o:s:' option; do
  case "$option" in
    a)
        annotation_gtf="$OPTARG"
        ;;
    e)
        ending="$OPTARG"
        ;;
    f)
        fastq_dir="$OPTARG"
        ;;
    g)
        genome_fasta="$OPTARG"
        ;;
    i)
        index=true
        ;;
    m)
        align_mode="$OPTARG"
        ;;
    o)
        output_star="$OPTARG"
        ;;  
    s)
        star_index="$OPTARG"
        ;;  
    h) 
        echo "$usage"
        exit 1
        ;;
   \?) 
        printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
done

# --- Check required options ---
if [[ -z "$star_index" || -z "$output_star" || -z "$align_mode" || -z "$genome_fasta" || -z "$annotation_gtf" || -z "$fastq_dir" || -z "$ending" ]]; then
  echo "Error: all -a -e -f -g -m -o -s options are required with arguments" >&2
  exit 1
fi

echo "All required options provided: Index directory : "$star_index", Output directory "$output_star", mode "$align_mode", genome FASTA "$genome_fasta", annotation "$annotation_gtf", FASTQ dir "$fastq_dir", file ending "$ending" "



if [ ! -d ${star_index} ]; then
    mkdir ${star_index}
fi

if [ ! -d ${output_star} ]; then
    mkdir ${output_star}
fi

if [[ "$index" == true ]]; then
    STAR --runThreadN 50 --runMode genomeGenerate --genomeDir $star_index --genomeFastaFiles ${genome_fasta}\
    --sjdbGTFfile $annotation_gtf --sjdbOverhang 49
    # --sjdbOverhang 49: maxreadlength - 1: this should actually be 37 for the leukemia samples, but for the others 30 should be fine
fi

shopt -s nullglob  # Prevents literal pattern if no match
files=($fastq_dir/*.fastq)

if [ ${#files[@]} -eq 0 ]; then
    cd $fastq_dir
    gunzip *.gz
    cd -
fi

# First try matching *R1* files
files=("${fastq_dir}"/*R1*.fastq)

if [ ${#files[@]} -eq 0 ]; then
    files=("${fastq_dir}"/*.fastq.1)
fi

if [ ${#files[@]} -eq 0 ]; then
    files=("${fastq_dir}"/*.fastq)
fi


for FQ in "${files[@]}"
do
    sample=$(basename "$FQ")       # remove path
    if [[ ${sample} == *"R1"* ]]; then
        sample=${sample%%R1*}          # remove R1 and everything after
    else
        sample=${sample%%.fastq*}
    fi
    echo ${sample}
    echo ${FQ}

    STAR\
    --runThreadN 16\
    --alignEndsType ${align_mode}\
    --outSAMstrandField intronMotif\
    --alignIntronMin 20\
    --alignIntronMax 1000000\
    --genomeDir $star_index\
    --readFilesIn ${FQ}\
    --twopassMode Basic\
    --seedSearchStartLmax 20\
    --outFilterMatchNminOverLread 0.9\
    --outSAMattributes All\
    --outSAMtype BAM SortedByCoordinate\
    --outFileNamePrefix "${output_star}"/"${sample}"${ending}

    samtools index --threads 32 "${output_star}"/"${sample}"${ending}Aligned.sortedByCoord.out.bam

    samtools idxstats "${output_star}"/"${sample}"${ending}Aligned.sortedByCoord.out.bam > \
    "${output_star}"/"${sample}"${ending}idxstats.out

    samtools stats "${output_star}"/"${sample}"${ending}Aligned.sortedByCoord.out.bam > \
    "${output_star}"/"${sample}"${ending}stats.out
done

# --alignMatesGapMax 20\# maximal genomic distance between mates, would like to set this to a small values as RPFs should fully overlap

# --seedSearchStartLmax 20\ reads will be split in pieces of 20 bp for MMP

# --peOverlapNbasesMin Minimum overlap of mates should be at least 30 bp, could also be set to 25 as this is the minimum
# minimum number of overlap bases to trigger mates merging and realignment

# --peOverlapMMp proportion of mismatches for overlapping sequence of the mates

# --alignEndsType Extend5pOfReads12: allow soft-clipping on 3' ends of both reads

# --outFilterMultimapNmax: defualt is 10, allow more in order to keep tRNA and rRNA alignments
#     --outFilterMultimapNmax 20\ makes files big and deduplication slow