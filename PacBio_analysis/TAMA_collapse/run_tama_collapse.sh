#!/bin/bash

#----- This script runs tama merge on 2 assemblies to merge them into one GTF file ----- #

eval "$(conda shell.bash hook)"
conda activate pacbio


reference_gtf="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"

bam_file_dir="/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment/HUVEC/pbmm2_align"

outdir_tama="/projects/splitorfs/work/PacBio/merged_bam_files/tama_collapse"

tama_tool_path="/home/ckalk/tools/tama"

tama_script_path="/home/ckalk/scripts/SplitORFs/PacBio_analysis/TAMA_collapse"

cell_type=HUVEC
if [[ ! -d "/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment/"${cell_type}"/pbmm2_align" ]]; then
    bash /home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/map_conditions/genome_mapping_cell_type.sh \
    -o "/projects/splitorfs/work/PacBio/merged_bam_files/genome_alignment" \
    -f "/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa" \
    -i "/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine" \
    -c "$cell_type"
fi

bash tama_collapse.sh -b "${bam_file_dir}"  -c "${cell_type}" -f "${genome_fasta}"\
 -o "${outdir_tama}" -r "${reference_gtf}" -t "${tama_tool_path}"


bash tama_merge_split_files.sh -t "${tama_tool_path}" -o "${outdir_tama}" \
 -c "${cell_type}" -s "${tama_script_path}"

 bash tama_merge_conditions.sh -t "${tama_tool_path}" -o "${outdir_tama}" \
 -c "${cell_type}" -s "${tama_script_path}"