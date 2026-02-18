#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

# =============================================================================
# Script Name: bowtie1_align.sh
# Description: This script aligns Ribo-seq reads end to end to the transcriptome
#              using bowtie1
# Usage:       bash bowtie1_align.sh args...
# Author:      Christina Kalk
# Date:        2025-10-10
# =============================================================================


bowtie1_base_name=$1
bowtie1_ref_fasta=$2
fastq_inpath=$3
bowtie1_out_dir=$4
aligned_name=$5
log_file=$6
script_dir=$7

if [ ! -d $bowtie1_out_dir ]; then
    mkdir $bowtie1_out_dir
fi

if [[ ${bowtie1_ref_fasta} != "no_index" ]]; then
    bowtie-build ${bowtie1_ref_fasta} ${bowtie1_base_name} --threads 32
    echo "done indexing"
fi

# First try matching *R1* files
shopt -s nullglob  # So non-matching globs result in empty arrays
files=("${fastq_inpath}"/*R1*.fastq.gz)

if [ ${#files[@]} -eq 0 ]; then
    files=("${fastq_inpath}"/*.fastq.1.gz)
fi

if [ ${#files[@]} -eq 0 ]; then
    files=("${fastq_inpath}"/*.fastq*)
fi

################################################################################
# Align files and filter alignments                                            #
################################################################################
echo "start alignment"
for FQ in "${files[@]}"
do
    sample=$(basename "$FQ")
    sample=${sample%%R1*}          # remove R1 and everything after

    outname="${bowtie1_out_dir}"/"${sample}"bowtie1_${aligned_name}
    echo ${sample}
    echo ${FQ}
    echo ${outname}

    # inspired from the ignolia paper
    bowtie \
        -v 2 \
        -m 200 \
        --norc \
        --best \
        -k 1 \
        -S \
        -q \
        --un ${outname}_unaligned.fastq \
        -x ${bowtie1_base_name} \
        ${FQ} \
        ${outname}_k1_R1.sam


    samtools view -@ 32 -bS ${outname}_k1_R1.sam \
     > ${outname}_k1_R1.bam

    samtools sort -@ 32 -o ${outname}_k1_R1_sorted.bam \
    ${outname}_k1_R1.bam

    samtools index -@ 32 ${outname}_k1_R1_sorted.bam

    samtools idxstats ${outname}_k1_R1_sorted.bam > \
    ${outname}_idxstats.out

    if [ ! -d $bowtie1_out_dir/filtered ]; then
        mkdir $bowtie1_out_dir/filtered
    fi
    filtered_name="${bowtie1_out_dir}"/filtered/"${sample}"bowtie1_${aligned_name}


    # remove secondary and supplementary alignments and unmapped
    samtools view -F 256 -F 2048 -F 0x4 -b ${outname}_k1_R1_sorted.bam > ${filtered_name}_k1_R1_sorted_filtered.bam
   

    samtools index -@ 32 ${filtered_name}_k1_R1_sorted_filtered.bam

    samtools idxstats ${filtered_name}_k1_R1_sorted_filtered.bam > \
    ${filtered_name}_idxstats.out

    samtools stats ${filtered_name}_k1_R1_sorted_filtered.bam > \
    ${filtered_name}_stats.out

    rm ${outname}_k1_R1.sam

done
# # report only 1 alignment (-k 1 is the default)
# # high penalty for mismatching bases in the seed
# # -m 200: allow up to 200 mapping positions



 python  "${script_dir}"/alignments/analyze_mappings/analyze_bowtie1_mappings.py \
 $log_file \
 ${bowtie1_out_dir}/summarized_bowtie_mapping_precents.csv


# ################################################################################
# # Filter for q10 alignments: filters rRNA and tRNA contmaination               #
# ################################################################################

# filter for quality of 10, this wil remove multi-mappers, hence contaminants
if [ ! -d $bowtie1_out_dir/filtered/q10 ]; then
    mkdir $bowtie1_out_dir/filtered/q10
fi

bam_files=("${bowtie1_out_dir}"/filtered/*.bam)
for bam in "${bam_files[@]}"
do
    sample=$(basename "$bam" .bam)

    # bowtie assigns 255 to uniquely mapping reads
    samtools view -q 255 -b $bam > \
     $bowtie1_out_dir/filtered/q10/${sample}_q10.bam

    samtools index $bowtie1_out_dir/filtered/q10/${sample}_q10.bam

    samtools idxstats $bowtie1_out_dir/filtered/q10/${sample}_q10.bam > \
    $bowtie1_out_dir/filtered/q10/${sample}_q10_idxstats.out

done


if [[ "$bowtie1_out_dir" != *Chae* ]]; then
    python  "${script_dir}"/alignments/summarize_bowtie2_alns_by_source.py \
    /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta \
    $bowtie1_out_dir/filtered

    python  "${script_dir}"/alignments/summarize_bowtie2_alns_by_source.py \
    /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta \
    $bowtie1_out_dir

    python  "${script_dir}"/alignments/summarize_bowtie2_alns_by_source.py \
    /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta \
    $bowtie1_out_dir/filtered/q10
 fi
