#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pacbio


bam_dir=$1
isoseq_outdir=$2
ref_fasta=$3
primer_file=$4

if [[ ! -d "$isoseq_outdir" ]]; then
    mkdir $isoseq_outdir
fi

if [[ ! -d "$isoseq_outdir"/refine ]]; then
    mkdir $isoseq_outdir/refine
fi

if [[ ! -d "$isoseq_outdir"/cluster ]]; then
    mkdir $isoseq_outdir/cluster
fi

if [[ ! -d "$isoseq_outdir"/mapped ]]; then
    mkdir $isoseq_outdir/mapped
fi

if [[ ! -d "$isoseq_outdir/collapsed" ]]; then
    mkdir $isoseq_outdir/collapsed
fi


shopt -s nullglob
bam_files=("${bam_dir}"/*bam)

ls "${bam_dir}"/HUVEC*bam > "${bam_dir}"/HUVEC_flnc.fofn

ls "${bam_dir}"/CM*bam > "${bam_dir}"/CM_flnc.fofn

for bam in "${bam_files[@]}"; 
do
    sample_name=$(basename $bam .primer_5p--primer_3p.bam)
    # isoseq refine --require-polya $bam primers.fasta $isoseq_outdir/refine/${sample_name}_refined.bam
    # I think require polyA will filter out reads wo polyA and anyway I do not have the primers file
    isoseq refine $bam $primer_file $isoseq_outdir/refine/${sample_name}_refined.bam --require-polya
done


#################################################################################
# ------------------ CLUSTERING                              ------------------ #
#################################################################################
export OMP_NUM_THREADS=32

# isoseq cluster2 "${bam_dir}"/HUVEC_flnc.fofn $isoseq_outdir/cluster/HUVEC_clustered.bam
# isoseq cluster2 "${bam_dir}"/CM_flnc.fofn $isoseq_outdir/cluster/CM_clustered.bam


#################################################################################
# ------------------ GENOME ALIGNEMENT PBMM2                 ------------------ #
#################################################################################
# pbmm2 align --preset ISOSEQ --sort $isoseq_outdir/cluster/HUVEC_clustered.bam $ref_fasta $isoseq_outdir/mapped/HUVEC_mapped.bam &
# pbmm2 align --preset ISOSEQ --sort $isoseq_outdir/cluster/CM_clustered.bam $ref_fasta $isoseq_outdir/mapped/CM_mapped.bam &
# wait


#################################################################################
# ------------------ ISOSEQ COLLAPSE                         ------------------ #
#################################################################################
# isoseq collapse $isoseq_outdir/mapped/HUVEC_mapped.bam $isoseq_outdir/collapsed/HUVEC_collapsed.gff &
# isoseq collapse $isoseq_outdir/mapped/CM_mapped.bam $isoseq_outdir/collapsed/CM_collapsed.gff &
# wait