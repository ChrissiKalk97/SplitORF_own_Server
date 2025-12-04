#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio

genome_fasta_file="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
ensembl_gtf_filtered="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"

pacbio_raw_datadir="/projects/splitorfs/work/own_data/PacBio_long_reads/run_23_06_25/objectstorage.uk-london-1.oraclecloud.com/Data-X208SC25032329-Z01-F001"
pacbio_merged_bamdir="/projects/splitorfs/work/PacBio/merged_bam_files"
raw_fastqc_dir=${pacbio_merged_bamdir}/fastqc
raw_longqc_dir=${pacbio_merged_bamdir}/longqc
lima_outdir=${pacbio_merged_bamdir}/lima
isoseq_outdir=${pacbio_merged_bamdir}/isoseq

primer_fasta=/projects/splitorfs/work/PacBio/merged_bam_files/lima/primer.fasta




#################################################################################
# ------------------ Pre-process PacBio long reads           ------------------ #
#################################################################################
if [[ ! -d "$pacbio_merged_bamdir" ]]; then
    mkdir $pacbio_merged_bamdir
fi

# remove the two folders with the failed samples
if [ -e /projects/splitorfs/work/own_data/PacBio_long_reads/run_23_06_25/objectstorage.uk-london-1.oraclecloud.com/Data-X208SC25032329-Z01-F001/HUVEC_DHYPO ]; then
    rm -r /projects/splitorfs/work/own_data/PacBio_long_reads/run_23_06_25/objectstorage.uk-london-1.oraclecloud.com/Data-X208SC25032329-Z01-F001/HUVEC_DHYPO
fi

if [ -e /projects/splitorfs/work/own_data/PacBio_long_reads/run_23_06_25/objectstorage.uk-london-1.oraclecloud.com/Data-X208SC25032329-Z01-F001/HUVEC_DMSO ]; then
    rm -r /projects/splitorfs/work/own_data/PacBio_long_reads/run_23_06_25/objectstorage.uk-london-1.oraclecloud.com/Data-X208SC25032329-Z01-F001/HUVEC_DMSO
fi


shopt -s nullglob
files=( "${pacbio_merged_bamdir}"/*merged.bam )

if (( ${#files[@]} == 0 )); then
    bash merge_pacbio_bams.sh "${pacbio_raw_datadir}" "${pacbio_merged_bamdir}"

    # FASTQC
    bash run_pacbio_fastqc.sh $pacbio_merged_bamdir $raw_fastqc_dir multiqc
fi

if [[ ! -e ""${pacbio_merged_bamdir}"/HUVEC_DHYPO_merged.bam" ]];then
    # copy the newly sequenced samples that do not need merging to the respective folder and name accordingly
    cp /projects/splitorfs/work/own_data/PacBio_long_reads/run_21_11_25/HUVEC_DHYPO.hifi_reads.bam \
    "${pacbio_merged_bamdir}"/HUVEC_DHYPO_merged.bam
    cp /projects/splitorfs/work/own_data/PacBio_long_reads/run_21_11_25/HUVEC_DHYPO.hifi_reads.bam.pbi \
    "${pacbio_merged_bamdir}"/HUVEC_DHYPO_merged.bam.pbi

    cp /projects/splitorfs/work/own_data/PacBio_long_reads/run_21_11_25/HUVEC_DMSO.hifi_reads.bam \
    "${pacbio_merged_bamdir}"/HUVEC_DMSO_merged.bam
    cp /projects/splitorfs/work/own_data/PacBio_long_reads/run_21_11_25/HUVEC_DMSO.hifi_reads.bam.pbi \
    "${pacbio_merged_bamdir}"/HUVEC_DMSO_merged.bam.pbi
fi

# LongQC
# bash run_longqc.sh ${pacbio_merged_bamdir} ${raw_longqc_dir}
# not running there are many errors

if [ ! -e $lima_outdir ];then
    # run lima for primer and polyA removal
    bash run_lima.sh $pacbio_merged_bamdir $lima_outdir $primer_fasta
fi

if [ ! -e $isoseq_outdir ];then
    # run isoseq pipeline steps
    bash run_isoseq.sh $lima_outdir $isoseq_outdir $genome_fasta_file $primer_fasta

    python helper_functions/plot_insert_lengths_after_refine.py \
    $isoseq_outdir/refine
fi

#################################################################################
# ------------------ Pre-process short reads                 ------------------ #
#################################################################################
# cd /home/ckalk/scripts/SplitORFs/short_RNA_seq

# raw_data_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData"
# merged_data_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/merged"
# raw_data_fastqc_dir_huvec="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/X208SC25032334-Z01-F001/01.RawData/fastqc"
# outidr_fastp="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_fastp"


# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/analyze_short_RNA_seq_cell_type.sh \
#  "${raw_data_dir_huvec}" \
#  "${merged_data_dir_huvec}" \
#  "${raw_data_fastqc_dir_huvec}" \
#  "${outidr_fastp}"

# raw_data_dir_cm="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/CM_short_reads/X208SC25032333-Z01-F002/01.RawData"
# merged_data_dir_cm="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/CM_short_reads/merged"
# raw_data_fastqc_dir_cm="/projects/splitorfs/work/own_data/Novogene/Michi_Vlado_run_1/CM_short_reads/X208SC25032333-Z01-F002/01.RawData/fastqc"
# outdir_fastp_cm="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/CM_fastp"


# bash /home/ckalk/scripts/SplitORFs/short_RNA_seq/analyze_short_RNA_seq_cell_type.sh \
#  "${raw_data_dir_cm}" \
#  "${merged_data_dir_cm}" \
#  "${raw_data_fastqc_dir_cm}" \
#  "${outdir_fastp_cm}"

# cd /home/ckalk/scripts/SplitORFs/PacBio_analysis


#################################################################################
# ------------------ Run Mandalorion                         ------------------ #
#################################################################################
# bash run_mandalorion_updated_parameters_correct_ref_one_cell_type_26_11_25.sh HUVEC
# bash run_mandalorion_updated_parameters_correct_ref_one_cell_type_26_11_25.sh CM




#################################################################################
# ------------------ Run Stringtie3                          ------------------ #
#################################################################################

# first need to align to the genome
# bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/map_conditions/map_conditions_to_assemblies.sh

# bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/run_stringtie3_correct_ref_06_10_25.sh


#################################################################################
# ------------------ TAMA merge final annotations            ------------------ #
#################################################################################

bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/compare_stringtie_mando/merge_stringtie_mando_correct_ref_06_10_25.sh