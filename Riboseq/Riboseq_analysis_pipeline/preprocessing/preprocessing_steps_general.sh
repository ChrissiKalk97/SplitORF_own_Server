#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq


# Help message:
usage="
Usage: ./Riboseq_pipeline.sh [-options] [arguments]

where:
-h			show this help
"
# available options for the programm
while getopts 'a:d:m:r:s:t:p:h' option; do
  case "$option" in
    a)
        adapter_sequence="$OPTARG"
        ;;
    d)
        data_dir="$OPTARG"
        ;;
    m) 
        max_length_given=true
        max_length="$OPTARG"
        ;;
    p) 
        fastp_dir="$OPTARG"
        ;;
    r)
        raw_fastqc_dir="$OPTARG"
        ;;
    s)
        script_dir="$OPTARG"
        ;;
    t)
        trim_tail_bool=true
        trim_tail="-$OPTARG"
        ;;
    h) 
        echo "$usage"
        exit 1
        ;;
   \?) 
        printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
done




################################################################################
# RUN FASTQ                                                                    #
################################################################################
# data_dir=$1
# raw_fastqc_dir=$2
# fastp_dir=$3
fastp_fastqc_dir=$raw_fastqc_dir/fastp
# adapter_sequence=$4

outdir_cutadapt=$(dirname "$raw_fastqc_dir")/cutadapt

if [[ ! "$max_length_given" == true ]]; then
    max_length=35
fi

if [[ ! -d  "${raw_fastqc_dir}" ]]; then
    mkdir "${raw_fastqc_dir}"
fi

if [[ ! -d  "${outdir_cutadapt}" ]]; then
    mkdir "${outdir_cutadapt}"
fi

if [[ ! -d  "${fastp_dir}" ]]; then
    mkdir "${fastp_dir}"
fi

if [[ ! -d  "${fastp_fastqc_dir}" ]]; then
    mkdir "${fastp_fastqc_dir}"
fi

shopt -s nullglob
samples=("${data_dir}"/*fastq)

################################################################################
# RUN FASTQC                                                                   #
################################################################################
for i in "${samples[@]}"
do
    fq=$i
    sample=$(basename "$i" .fastq)
    if [[ ! -e  "${raw_fastqc_dir}/${sample}_fastqc.html" ]]; then
    (
        fastqc \
        -o "${raw_fastqc_dir}"/ \
        -t 32\
        ${fq} 
        )&
    fi
done

wait

if [[ ! -e  "${raw_fastqc_dir}/fastqc_multiqc.html" ]]; then
    bash /home/ckalk/scripts/SplitORFs/Riboseq/preprocess_data/fastqc_multiqc_for_all.sh \
        ${data_dir} \
        "${raw_fastqc_dir}" \
        fastqc_multiqc \
        FASTQ
fi



for fq in "${samples[@]}"
do
    sample=$(basename "$fq" .fastq)
    if [[ ! -e  "${fastp_dir}/${sample}.cutadapt.fastq" ]]; then
        if [[ $trim_tail_bool == true ]]; then
            cutadapt \
                -a ${adapter_sequence} \
                -o "${outdir_cutadapt}"/${sample}.cutadapt.fastq \
                --cores 16 \
                --minimum-length 25 \
                --max-n 0.1 \
                --max-expected-errors 1 \
                -u ${trim_tail} \
                ${fq}
        else
            cutadapt \
                -a ${adapter_sequence} \
                -o "${outdir_cutadapt}"/${sample}.cutadapt.fastq \
                --cores 16 \
                --minimum-length 25 \
                --max-n 0.1 \
                --max-expected-errors 1 \
                ${fq}
        fi
    fi
    if [[ ! -e  "${fastp_dir}/${sample}_fastp.fastq" ]]; then
        # echo "${fastp_dir}"/${sample}_fastp.fastq
        fastp \
            -i "${outdir_cutadapt}"/${sample}.cutadapt.fastq \
            -o "${fastp_dir}"/${sample}_fastp.fastq \
            --json "${fastp_dir}"/${sample}.fastp.json \
            --thread 32 \
            --length_required 25 \
            --cut_front \
            --cut_front_window_size 1\
            --cut_mean_quality 15 \
            --cut_tail \
            --cut_tail_window_size 1\
            --length_limit ${max_length}
    fi
done





# ################################################################################
# # RUN FASTP                                                                    #
# ################################################################################
# samples=("${outdir_cutadapt}"/*fastq)

# for i in "${samples[@]}"
# do
#     sample=$(basename "$i" .cutadapt.fastq)
#     fq=$i

#     echo $fq

#     if [[ ! -e  "${fastp_dir}/${sample}_fastp.fastq" ]]; then

#         # echo "${fastp_dir}"/${sample}_fastp.fastq
#         fastp \
#             -i ${fq} \
#             -o "${fastp_dir}"/${sample}_fastp.fastq \
#             --json "${fastp_dir}"/${sample}.fastp.json \
#             --thread 32 \
#             --length_required 25 \
#             --cut_front \
#             --cut_front_window_size 1\
#             --cut_mean_quality 15 \
#             --cut_tail \
#             --cut_tail_window_size 1\
#             --length_limit ${max_length}
#     fi
# done



# ################################################################################
# # RUN FASTQC AFTER FASTP                                                       #
# ################################################################################
bash /home/ckalk/scripts/SplitORFs/Riboseq/preprocess_data/fastqc_multiqc_for_all.sh \
    "${fastp_dir}" \
    "${fastp_fastqc_dir}" \
    fastp_multiqc \
    FASTQ