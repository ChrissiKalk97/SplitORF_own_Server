#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

UMI_Indir=$1
UMI_dedup_outdir=$2
align_type=$3

shopt -s nullglob
if [[ "$align_type" == "transcriptomic" ]]; then
    echo $align_type
    files=("${UMI_Indir}"/*filtered*.bam)
else
    echo $align_type
    files=("${UMI_Indir}"/*.bam)
fi



for BAM in "${files[@]}"
do
    (
        filename=$(basename "$BAM")
        # are the bamfiles genome or transcriptome aligned
        # remove respective suffix to get the sample name
        if [[ "$filename" == *.cutadapt_umi_fastp.only_R1_Aligned.sortedByCoord.out.bam ]]; then
            sample=${filename%.cutadapt_umi_fastp.only_R1_Aligned.sortedByCoord.out.bam}
        else
            if [[ "$filename" == *"q10"* ]]; then
                sample=${filename%.cutadapt_umi_fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered_q10.bam}  
            else
                sample=${filename%.cutadapt_umi_fastp.bowtie2_concat_transcriptome_k1_R1_sorted_filtered.bam}  
            fi
        fi
        
        echo $sample
        umi_tools dedup \
        --stdin=${BAM} \
        --log="${UMI_dedup_outdir}"/${sample}_LOGFILE \
        --output-stats="${UMI_dedup_outdir}"/${sample}_outstats \
            > "${UMI_dedup_outdir}"/${sample}_dedup.bam

        samtools index "${UMI_dedup_outdir}"/${sample}_dedup.bam

        samtools flagstat "${UMI_dedup_outdir}"/${sample}_dedup.bam > \
        "${UMI_dedup_outdir}"/"${sample}"_dedup_flagstat.out

        samtools idxstats  "${UMI_dedup_outdir}"/${sample}_dedup.bam > \
        "${UMI_dedup_outdir}"/"${sample}".dedup_idxstats.out
    ) &
done

wait

python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/summarize_bowtie2_alns_by_source.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta \
 $UMI_dedup_outdir