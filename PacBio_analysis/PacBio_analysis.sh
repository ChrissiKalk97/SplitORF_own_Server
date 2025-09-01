#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

genome_fasta_file="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
ensembl_gtf_filtered="/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"

pacbio_raw_datadir="/projects/splitorfs/work/own_data/PacBio_long_reads/run_23_06_25/objectstorage.uk-london-1.oraclecloud.com/Data-X208SC25032329-Z01-F001"
pacbio_merged_bamdir="/projects/splitorfs/work/PacBio/merged_bam_files"
raw_fastqc_dir=${pacbio_merged_bamdir}/fastqc
raw_longqc_dir=${pacbio_merged_bamdir}/longqc
lima_outdir=${pacbio_merged_bamdir}/lima
isoseq_outdir=${pacbio_merged_bamdir}/isoseq

primer_fasta=/projects/splitorfs/work/PacBio/merged_bam_files/lima/primer.fasta

if [[ ! -d "$pacbio_merged_bamdir" ]]; then
    mkdir $pacbio_merged_bamdir
fi

# bash merge_pacbio_bams.sh "${pacbio_raw_datadir}" "${pacbio_merged_bamdir}"


# FASTQC
# bash run_pacbio_fastqc.sh $pacbio_merged_bamdir $raw_fastqc_dir multiqc


# LongQC
# bash run_longqc.sh ${pacbio_merged_bamdir} ${raw_longqc_dir}
# not running there are many errors

# run lima for primer and polyA removal
# bash run_lima.sh $pacbio_merged_bamdir $lima_outdir $primer_fasta

# run isoseq pipeline steps
# bash run_isoseq.sh $lima_outdir $isoseq_outdir $genome_fasta_file $primer_fasta

python helper_functions/plot_insert_lengths_after_refine.py \
 $isoseq_outdir/refine
