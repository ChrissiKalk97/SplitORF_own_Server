#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

UMI_Indir=$1
UMI_dedup_outdir=$2


for BAM in "${UMI_Indir}"/*.bam

do
    sample=$(basename ${BAM} .cutadapt_umi_fastp.only_R1_Aligned.sortedByCoord.out.bam)
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
done


python /home/ckalk/scripts/SplitORFs/Riboseq/Michi_Vlado_round_1/alignments/summarize_bowtie2_alns_by_source.py \
 /projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mtDNA/Homo_sapiens.GRCh38.dna.chromosome.MT.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/rRNA/redownload/rRNA_ref_NCBI_Ens.fasta,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/tRNA/hg38-tRNAs/hg38-tRNAs.fa,/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ens_Gencode_lncRNA_ncRNA.fasta \
 $UMI_dedup_outdir