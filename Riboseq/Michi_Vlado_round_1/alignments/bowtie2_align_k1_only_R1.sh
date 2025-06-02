#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq


Bowtie2_base_name=$1
Bowtie2_ref_fasta=$2
FASTQ_Inpath=$3
Bowtie2_out_dir=$4
aligned_name=$5
log_file=$6

if [[ ${Bowtie2_ref_fasta} != "no_index" ]]; then
    bowtie2-build ${Bowtie2_ref_fasta} ${Bowtie2_base_name} --threads 32
fi

# First try matching *R1* files
shopt -s nullglob  # So non-matching globs result in empty arrays
files=("${FASTQ_Inpath}"/*R1*.fastq.gz)

if [ ${#files[@]} -eq 0 ]; then
    files=("${FASTQ_Inpath}"/*.fastq.1.gz)
fi


for FQ in "${files[@]}"
do
    sample=$(basename "$FQ")
    sample=${sample%%R1*}          # remove R1 and everything after

    outname="${Bowtie2_out_dir}"/"${sample}"bowtie2_${aligned_name}
    echo ${sample}
    echo ${FQ}
    echo ${outname}



    bowtie2 \
    -p 32 \
    --np 3 \
    --very-sensitive-local \
    --ignore-quals \
    --score-min L,0,1.8 \
    -S ${outname}_k1_R1.sam \
    --un-conc-gz ${outname}_unaligned.fastq.gz \
    -x ${Bowtie2_base_name} \
    -q \
    -U ${FQ}

    samtools view -@ 32 -bS ${outname}_k1_R1.sam \
     > ${outname}_k1_R1.bam

    samtools sort -@ 32 -o ${outname}_k1_R1_sorted.bam \
    ${outname}_k1_R1.bam


    if [ ! -d $Bowtie2_out_dir/filtered ]; then
        mkdir $Bowtie2_out_dir/filtered
    fi
    filtered_name="${Bowtie2_out_dir}"/filtered/"${sample}"bowtie2_${aligned_name}

#     # remove secondary and supplementary alignments
    samtools view -F 256 -F 2048 -F 0x4 -b ${outname}_k1_R1_sorted.bam > \
     ${filtered_name}_k1_R1_sorted_filtered.bam
   

    samtools index --threads 32 ${filtered_name}_k1_R1_sorted_filtered.bam

    samtools idxstats ${filtered_name}_k1_R1_sorted_filtered.bam > \
    ${filtered_name}_idxstats.out

    samtools stats ${filtered_name}_k1_R1_sorted_filtered.bam > \
    ${filtered_name}_stats.out



    rm ${outname}_k1_R1.sam
done
# # report only 1 alignment (-k 1 is the default)
# # high penalty for mismatching bases in the seed
# # -N 1 \ allow for one mismatch in the seed alignment, because of the polyA
# # -L 20 \   # 20 is actually the default in local mode


python  /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/summarize_bowtie2_alns_by_source.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta \
 $Bowtie2_out_dir/filtered

 python  /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/preprocessing/analyze_mappings/analyze_bowtie2_mappings.py \
 $log_file \
 ${Bowtie2_out_dir}/summarized_bowtie_mapping_precents.csv

# filter for quality of 10, this wil remove multi-mappers, hence contaminants
if [ ! -d $Bowtie2_out_dir/filtered/q10 ]; then
    mkdir $Bowtie2_out_dir/filtered/q10
fi

bam_files=("${Bowtie2_out_dir}"/filtered/*.bam)
for bam in "${bam_files[@]}"
do
    sample=$(basename "$bam" .bam)

    samtools view -q 10 -b $bam > \
     $Bowtie2_out_dir/filtered/q10/${sample}_q10.bam

    samtools index $Bowtie2_out_dir/filtered/q10/${sample}_q10.bam

    samtools idxstats $Bowtie2_out_dir/filtered/q10/${sample}_q10.bam > \
    $Bowtie2_out_dir/filtered/q10/${sample}_q10_idxstats.out

done

python  /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/summarize_bowtie2_alns_by_source.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta \
 $Bowtie2_out_dir/filtered/q10