#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

################################################################################
# RUN FASTQ                                                                    #
################################################################################
data_dir=$1
fastqc_dir=$2
fastp_dir=$3
fastp_fastqc_dir=$fastqc_dir/fastp

# set maxlength to the 4th argument if it is set, else to 35
max_length=${4:-35}

if [[ -d  ${fastqc_dir} ]]; then
    mkdir $fastqc_dir
fi

if [[ -d  ${fastp_dir} ]]; then
    mkdir $fastp_dir
fi

if [[ -d  ${fastp_fastqc_dir} ]]; then
    mkdir $fastp_fastqc_dir
fi

shopt -s nullglob
samples=("${data_dir}"/*fastq)

################################################################################
# RUN FASTQC                                                                   #
################################################################################
# for i in "${samples[@]}"
# do
#     fq=$i
#     sample=$(basename "$i" .fastq)
#     if [ ! -e  "${fastqc_dir}/${sample}_fastqc.html" ]; then
#     (
#         fastqc \
#         -o ${fastqc_dir}/ \
#         -t 32\
#         ${fq} 
#         )&
#     fi
# done

# wait

# bash /home/ckalk/scripts/SplitORFs/Riboseq/preprocess_data/fastqc_multiqc_for_all.sh \
#     ${data_dir} \
#     ${fastqc_dir} \
#     fastqc_multiqc \
#     FASTQ


# ################################################################################
# # RUN FASTP                                                                    #
# ################################################################################
for i in "${samples[@]}"
do
    sample=$(basename "$i" .fastq)
    fq=$i

    echo $fq

    if [ ! -e  "${fastp_dir}/${sample}_fastp.fastq" ]; then

        echo ${fastp_dir}/${sample}_fastp.fastq
        fastp \
            -i ${fq} \
            -o ${fastp_dir}/${sample}_fastp.fastq \
            --json ${fastp_dir}/${sample}.fastp.json \
            --thread 32 \
            --length_required 25 \
            --length_limit ${max_length}
    fi
done



# ################################################################################
# # RUN FASTQC AFTER FASTP                                                       #
# ################################################################################
bash /home/ckalk/scripts/SplitORFs/Riboseq/preprocess_data/fastqc_multiqc_for_all.sh \
    ${fastp_dir} \
    ${fastp_fastqc_dir} \
    fastp_multiqc \
    FASTQ