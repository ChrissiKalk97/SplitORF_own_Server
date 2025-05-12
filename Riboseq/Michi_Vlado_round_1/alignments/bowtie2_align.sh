#!/bin/bash

Bowtie2_base_name=$1
Bowtie_ref_fasta=$2
FASTQ_Inpath=$3
Bowtie_out_dir=$4
aligned_name=$5

if [[ ${Bowtie_ref_fasta} != "no_index" ]]; then
    bowtie2-build ${Bowtie_ref_fasta} ${Bowtie2_base_name} --threads 32
fi

# First try matching *R1* files
shopt -s nullglob  # So non-matching globs result in empty arrays
files=("${FASTQ_Inpath}"/*R1*.fastq.gz)

if [ ${#files[@]} -eq 0 ]; then
    files=("${FASTQ_Inpath}"/*.fastq.1.gz)
fi


for FQ in "${files[@]}"
do
    sample=$(basename "$FQ")       # remove path
    if [[ ${sample} == *"R1"* ]]; then
        sample=${sample%%R1*}          # remove R1 and everything after
        FQ2=${FQ/R1/R2} # substitute R1 with R2 in the whole file path
    else
        sample=${sample%%.fastq*} 
        FQ2=${FQ/.1/.2} # substitute R1 with R2 in the whole file path
    fi
    echo ${sample}
    echo ${FQ}
    echo ${FQ2}
    bowtie2 \
    -p 32 \
    -k 100 \
    --np 3 \
    --very-sensitive-local \
    --no-discordant \
    --ignore-quals \
    --score-min L,0,1.8 \
    -S "${Bowtie_out_dir}"/"${sample}"bowtie2_${aligned_name}.sam \
    --un-conc-gz "${Bowtie_out_dir}"/"${sample}"bowtie2_${aligned_name}_unaligned.fastq.gz \
    -x ${Bowtie2_base_name} \
    -q \
    -1 ${FQ} \
    -2 ${FQ2}
done
# report up to 100 alignments
# high penalty for mismatching bases in the seed
# -N 1 \ allow for one mismatch in the seed alignment, because of the polyA
# -L 20 \   # 20 is actually the default in local mode