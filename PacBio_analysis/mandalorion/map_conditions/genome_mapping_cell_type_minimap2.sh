#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio


# out_path="/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment"
# genome_fasta_file="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

# isoseq_reads_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"



# Define the usage function
usage() {
  echo "Usage: $0 -i <isoseq_reads_dir> -f <genome_fasta> -g <gtf> -o <out_path>  -c <cell_type> [-h for help]"
}

# Process options with silent error mode
while getopts "i:f:g:o:c:h" opt; do
  case $opt in
    c)
      cell_type="$OPTARG"
      ;;
    f)
      genome_fasta_file="$OPTARG"
      ;;
    g)
      gtf="$OPTARG"
      ;;
    i)
      isoseq_reads_dir="$OPTARG"
      ;;
    o)
      out_path="$OPTARG"
      ;;
    h)
      usage
      exit 0
      ;;
    :)
      echo "Error: Option -$OPTARG requires an argument."
      usage
      exit 1
      ;;
    \?)
      echo "Error: Invalid option -$OPTARG"
      usage
      exit 1
      ;;
  esac
done







if [[ ! -d "$out_path" ]]; then
    mkdir "$out_path"
fi

if [[ ! -d "$out_path/"${cell_type}"" ]]; then
    mkdir "$out_path/"${cell_type}""
fi


if [[ ! -d "$out_path/"${cell_type}"/minimap2_align" ]]; then
    mkdir "$out_path/"${cell_type}"/minimap2_align"
fi

mkdir -p "$out_path/"${cell_type}"/minimap2_align/genome_index"
#################################################################################
# ------------------ ALIGN TO GENOME                         ------------------ #
#################################################################################
# might not want to index, as the parameters need to be adjusted for the datatype
# and I simply could not find the informaiton about how to do this 
# for my parameters

if [[ ! -e "$out_path/"${cell_type}/minimap2_align/genome_index/"$name_gtf"_junctions.bed ]]; then
    name_gtf=$(basename "$gtf" .gtf)
    paftools.js gff2bed "$gtf" > "$out_path/"${cell_type}/minimap2_align/genome_index/"$name_gtf"_junctions.bed
fi



for fastq in "${isoseq_reads_dir}"/"${cell_type}"*fastq; do
    sample=$(basename "$fastq")
    sample="${sample%_merged_lima_refined.fastq}"
    if [ ! -e "${out_path}"/"${cell_type}"/minimap2_align/"${sample}"_minimap2_aligned_genome_sorted.bam ]; then
        echo "${out_path}"/"${cell_type}"/minimap2_align/"${sample}"_minimap2_aligned_genome.bam
        minimap2 -ax splice:hq -uf --junc-bed="$out_path/${cell_type}/minimap2_align/genome_index/"$name_gtf"_junctions.bed" -t 40\
        "$genome_fasta_file" \
        "$fastq" | samtools sort -o \
        "${out_path}"/"${cell_type}"/minimap2_align/"${sample}"_minimap2_aligned_genome.bam
    fi
 done

shopt -s nullglob
bam_files=("${out_path}/${cell_type}/minimap2_align/"*aligned_genome.bam)

if (( ${#bam_files[@]} > 0 )); then
  for bam in "${out_path}"/"${cell_type}"/minimap2_align/*aligned_genome.bam; do
      if [[ ! -e "$(dirname "$bam")/$(basename "$bam" .bam)_sorted.bam" ]]; then
        samtools sort -@ 30 -o $(dirname "$bam")/$(basename "$bam" .bam)_sorted.bam "$bam"
        samtools index -@ 30 $(dirname "$bam")/$(basename "$bam" .bam)_sorted.bam
        samtools idxstats -@ 30 $(dirname "$bam")/$(basename "$bam" .bam)_sorted.bam > \
        $(dirname "$bam")/$(basename "$bam" .bam)_idxstats.out
        samtools flagstat -@ 30 $(dirname "$bam")/$(basename "$bam" .bam)_sorted.bam > \
        $(dirname "$bam")/$(basename "$bam" .bam)_flagstat.out
        samtools stats -@ 30 $(dirname "$bam")/$(basename "$bam" .bam)_sorted.bam > \
        $(dirname "$bam")/$(basename "$bam" .bam)_stats.out
        samtools view -@ 30 -q 20 -F 0x904 -b $(dirname "$bam")/$(basename "$bam" .bam)_sorted.bam | \
        samtools sort -@ 30 \
        -o $(dirname "$bam")/$(basename "$bam" .bam)_filtered.bam
      fi
  done
fi

