#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

################################################################################
# RUN FASTQ                                                                    #
################################################################################
OUTDIR="/projects/splitorfs/work/Riboseq/data/fastqc"
heart_dir="/projects/splitorfs/work/Riboseq/data/heart_iPSC/fastq"
leukemia_dir="/projects/splitorfs/work/Riboseq/data/leukemia/fastq"
endothel_dir="/projects/splitorfs/work/Riboseq/data/endothel_Siragusa/fastq_from_bam"

samples=($heart_dir/"ERR3367797.fastq" $heart_dir/"ERR3367798.fastq")

for i in {1..6}; do
  samples+=($endothel_dir"/OHMX20220060_00$i".fastq)
done 

for i in {8..9}; do
  samples+=($leukemia_dir"/SRR1129460$i".fastq)
done

for i in {10..11}; do
  samples+=($leukemia_dir"/SRR112946$i".fastq)
done

################################################################################
# RUN FASTQC                                                                   #
################################################################################
# for i in "${samples[@]}"
# do
# FQ=$i
# fastqc \
#     -o ${OUTDIR}/ \
#     -t 32\
#     ${FQ}
# done


################################################################################
# RUN FASTP                                                                    #
################################################################################
OUTDIR="/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples"

for i in "${samples[@]}"
do
SAMPLE=$(basename "$i" .fastq)
FQ=$i
fastp \
    -i ${FQ} \
    -o ${OUTDIR}/${SAMPLE}_fastp.fastq \
    --json ${OUTDIR}/${SAMPLE}.fastp.json \
    --thread 32 \
    --length_required 25 \
    --length_limit 35
done



################################################################################
# RUN FASTQC AFTER FASTP                                                       #
################################################################################
OUTDIR="/projects/splitorfs/work/Riboseq/data/fastqc/fastp"
fastp_dir="/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples"

samples=()

for fastq in $fastp_dir/*.fastq; do
  samples+=($fastq)
done 

for i in "${samples[@]}"
do
FQ=$i
fastqc \
    -o ${OUTDIR}/ \
    -t 32\
    ${FQ}
done
