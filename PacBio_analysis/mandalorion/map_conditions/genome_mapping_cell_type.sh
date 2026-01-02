#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio


# out_path="/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment"
# genome_fasta_file="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

# isoseq_reads_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"



# Define the usage function
usage() {
  echo "Usage: $0 -i <isoseq_reads_dir> -g <genome_fasta> -o <out_path>  -c <cell_type> [-h for help]"
}

# Process options with silent error mode
while getopts "i:f:o:c:h" opt; do
  case $opt in
    i)
      isoseq_reads_dir="$OPTARG"
      ;;
    f)
      genome_fasta_file="$OPTARG"
      ;;
    o)
      out_path="$OPTARG"
      ;;
    c)
      cell_type="$OPTARG"
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


if [[ ! -d "$out_path/"${cell_type}"/pbmm2_align" ]]; then
    mkdir "$out_path/"${cell_type}"/pbmm2_align"
fi


if [[ ! -d "$out_path/"${cell_type}"/pbmm2_align/genome_index" ]]; then
    mkdir "$out_path/"${cell_type}"/pbmm2_align/genome_index"
fi


#################################################################################
# ------------------ ALIGN TO GENOME                         ------------------ #
#################################################################################
if [ ! -e "$out_path/"${cell_type}"/pbmm2_align/genome_index/genome_index.mmi" ]; then 
    pbmm2 index ${genome_fasta_file} $out_path/"${cell_type}"/pbmm2_align/genome_index/genome_index.mmi --preset ISOSEQ
fi


for bam in "${isoseq_reads_dir}"/"${cell_type}"*bam; do
    sample=$(basename $bam)
    sample="${sample%_merged_lima_refined.bam}"
    if [ ! -e $out_path/"${cell_type}"/pbmm2_align/${sample}_pbmm2_aligned_genome.bam ]; then
        echo $out_path/"${cell_type}"/pbmm2_align/${sample}_pbmm2_aligned_genome.bam
        pbmm2 align $out_path/"${cell_type}"/pbmm2_align/genome_index/genome_index.mmi \
        $bam \
        $out_path/"${cell_type}"/pbmm2_align/${sample}_pbmm2_aligned_genome.bam \
        --preset ISOSEQ
    fi
 done



for bam in $out_path/"${cell_type}"/pbmm2_align/*.bam; do
    if [ ! -e $(dirname $bam)/$(basename $bam .bam)_filtered.bam ]; then
      samtools sort -@ 30 -o $(dirname $bam)/$(basename $bam .bam)_sorted.bam $bam
      samtools index -@ 30 $(dirname $bam)/$(basename $bam .bam)_sorted.bam
      samtools idxstats -@ 30 $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
      $(dirname $bam)/$(basename $bam .bam)_idxstats.out
      samtools flagstat -@ 30 $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
      $(dirname $bam)/$(basename $bam .bam)_flagstat.out
      samtools stats -@ 30 $(dirname $bam)/$(basename $bam .bam)_sorted.bam > \
      $(dirname $bam)/$(basename $bam .bam)_stats.out
      samtools view -@ 30 -q 20 -F 0x904 -b $(dirname $bam)/$(basename $bam .bam)_sorted.bam | \
      samtools sort -@ 30 \
      -o $(dirname $bam)/$(basename $bam .bam)_filtered.bam
    fi
done

