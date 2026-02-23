#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

umi_indir=$1
umi_dedup_outdir=$2
align_type=$3
module_dir=$4


if [ ! -d ${umi_dedup_outdir} ]; then
    mkdir ${umi_dedup_outdir}
fi

shopt -s nullglob
if [[ "$align_type" == "transcriptomic" ]]; then
    echo $align_type
    files=("${umi_indir}"/*filtered*.bam)
    parallel=true
else
    echo $align_type
    files=("${umi_indir}"/*.bam)
    parallel=false
fi


if [[ "${parallel}" == true ]]; then
    for bam in "${files[@]}"
    do
        (
            filename=$(basename "$bam")
            # are the bamfiles genome or transcriptome aligned
            # remove respective suffix to get the sample name
            if [[ "$filename" == *"q10"*"bowtie2"* ]]; then
                sample=${filename%.cutadapt_umi_fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10.bam}  
            elif [[ "$filename" == *"bowtie2"* ]]; then
                sample=${filename%.cutadapt_umi_fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered.bam} 
            elif [[ "$filename" == *"q10"*"bowtie1"* ]]; then
                sample=${filename%.cutadapt_umi_fastp.bowtie1_concat_transcriptome_k1_R1_sorted_filtered_q10.bam}  
            else
                sample=${filename%.cutadapt_umi_fastp.bowtie1_concat_transcriptome_k1_R1_sorted_filtered.bam}  
            fi

            
            echo $sample
            umi_tools dedup \
            --stdin=${bam} \
            --log="${umi_dedup_outdir}"/${sample}_LOGFILE \
            --output-stats="${umi_dedup_outdir}"/${sample}_outstats \
                > "${umi_dedup_outdir}"/${sample}_dedup.bam

            samtools index "${umi_dedup_outdir}"/${sample}_dedup.bam

            samtools flagstat "${umi_dedup_outdir}"/${sample}_dedup.bam > \
            "${umi_dedup_outdir}"/"${sample}"_dedup_flagstat.out

            samtools idxstats  "${umi_dedup_outdir}"/${sample}_dedup.bam > \
            "${umi_dedup_outdir}"/"${sample}".dedup_idxstats.out
        ) &
    done

    wait
else
    for bam in "${files[@]}"
    do

            filename=$(basename "$bam")
            # are the bamfiles genome or transcriptome aligned
            # remove respective suffix to get the sample name
            if [[ "$filename" == *.cutadapt_umi_fastp.only_R1_Aligned.sortedByCoord.out.bam ]]; then
                sample=${filename%.cutadapt_umi_fastp.only_R1_Aligned.sortedByCoord.out.bam}
            elif [[ "$filename" == *_Aligned.sortedByCoord.out.bam ]]; then
                sample=${filename%_Aligned.sortedByCoord.out.bam}
            fi
            
            echo $sample
            umi_tools dedup \
            --stdin=${bam} \
            --log="${umi_dedup_outdir}"/${sample}_LOGFILE \
            --output-stats="${umi_dedup_outdir}"/${sample}_outstats \
                > "${umi_dedup_outdir}"/${sample}_dedup.bam

            samtools index "${umi_dedup_outdir}"/${sample}_dedup.bam

            samtools flagstat "${umi_dedup_outdir}"/${sample}_dedup.bam > \
            "${umi_dedup_outdir}"/"${sample}"_dedup_flagstat.out

            samtools idxstats  "${umi_dedup_outdir}"/${sample}_dedup.bam > \
            "${umi_dedup_outdir}"/"${sample}".dedup_idxstats.out
    done
fi

if [[ "$align_type" == "transcriptomic" ]]; then
    python "${module_dir}"/alignments/summarize_bowtie2_alns_by_source.py \
    /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta \
    $umi_dedup_outdir
fi