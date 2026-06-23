#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate isoquant

genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
ensembl_gtf_filtered="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"

output_dir="/projects/splitorfs/work/PacBio/merged_bam_files/IsoQuant"

bam_dir="/projects/splitorfs/work/PacBio/merged_bam_files/isoseq/refine"
script_dir="/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3"

stringtie3_dir_raw="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/IsoQuant_raw_June_2026"
sqanti_dir="$output_dir"/SQANTI3
sqanti_script_dir="/home/ckalk/scripts/SplitORFs/PacBio_analysis/mandalorion/assess_mando_sqanti3"


mkdir -p "$output_dir"


# assumpton: Preprocessing steps from 
# /home/ckalk/scripts/SplitORFs/PacBio_analysis/PacBio_analysis_29_05_26.sh
# are alaready run!

# might want to add this script as a subsript to the PacBio analysis



#################################################################################
# ------------------ Run IsoQuant                            ------------------ #
#################################################################################

# isoquant -d pacbio_ccs \
#  --polya_trimmed stranded \
#  --sqanti_output \
#  --yaml /home/ckalk/scripts/SplitORFs/PacBio_analysis/IsoQuant/long_reads.yaml  \
#  --complete_genedb --genedb "${ensembl_gtf_filtered}" \
#   --reference "$genome_fasta" --output "$output_dir"


 # --fl_data \
 # this option means that both ends of the reads are reliable, which is not the case


#################################################################################
# ------------------ Run SQANTIQC                            ------------------ #
#################################################################################
if [ ! -d "${stringtie3_dir_raw}" ]; then
    mkdir "${stringtie3_dir_raw}"
fi

if [ ! -d "${stringtie3_dir_raw}"/kallisto ]; then
    mkdir "${stringtie3_dir_raw}"/kallisto
fi

if [ ! -d "${stringtie3_dir_raw}"/kallisto/index ]; then
    mkdir "${stringtie3_dir_raw}"/kallisto/index
fi

for cell_type in CM HUVEC; do

    echo $cell_type
    short_read_dir="/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/${cell_type}_fastp"

    # if [ ! -d "${stringtie3_dir_raw}"/kallisto/"${cell_type}"_quant ]; then
    #     mkdir "${stringtie3_dir_raw}"/kallisto/"${cell_type}"_quant

    #     echo "$output_dir"/"${cell_type}"/"${cell_type}".transcript_models.gtf

    #     bash ${script_dir}/kallisto/kallisto_index.sh \
    #     "$output_dir"/"${cell_type}"/"${cell_type}".transcript_models.gtf \
    #     ${genome_fasta} \
    #     "${stringtie3_dir_raw}"/kallisto/"${cell_type}"_isoquant_assembly_transcriptome.fa \
    #     "${stringtie3_dir_raw}"/kallisto/index/"${cell_type}"

    #     bash ${script_dir}/kallisto/kallisto_quantification.sh \
    #     "${stringtie3_dir_raw}"/kallisto/index/"${cell_type}".idx \
    #     "${short_read_dir}" \
    #     "${stringtie3_dir_raw}"/kallisto/"${cell_type}"_quant
    # fi


    if [ ! -d "${sqanti_dir}" ]; then
        mkdir "${sqanti_dir}"
    fi


    if [ ! -d "${sqanti_dir}"/SQANTI3_QC ]; then
        mkdir "${sqanti_dir}"/SQANTI3_QC
    fi


    # if [ ! -d "${sqanti_dir}"/SQANTI3_QC/"${cell_type}" ]; then
        
        
        #mkdir -p "${sqanti_dir}"/SQANTI3_QC/"${cell_type}"
        # bash ${script_dir}/sqanti3/sqanti3_qc_mando_huvec.sh \
        # /home/ckalk/tools/sqanti3.6 \
        # "$output_dir"/"${cell_type}"/"${cell_type}".transcript_models.gtf \
        # ${ensembl_gtf_filtered} \
        # ${genome_fasta} \
        # "${sqanti_dir}"/SQANTI3_QC/"${cell_type}" \
        # ${script_dir}/sqanti3/"${cell_type}"_short_reads.txt \
        # "${stringtie3_dir_raw}"/kallisto/"${cell_type}"_quant \
        # "$output_dir"/"$cell_type"/"$cell_type"_isoquant_fl_for_sqanti.tsv

    bash "${script_dir}"/run_sqanti3_on_mando_one_cell_type_20_10_25.sh \
    ${cell_type} \
    "$genome_fasta" \
    "$ensembl_gtf_filtered" \
    "$output_dir" \
    "${stringtie3_dir_raw}" \
    "$short_read_dir" \
    "${sqanti_script_dir}" \
    "$bam_dir" \
    "isoquant"
    # fi
done