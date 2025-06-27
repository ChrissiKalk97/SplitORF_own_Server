#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate translationAI

cd /projects/splitorfs/work/LLMs/TranslationAI/Output

echo "Reformat fasta header"

python /home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/reformat_fasta_header.py\
 /projects/splitorfs/work/LLMs/TranslationAI/Input/RI_transcripts_110_TranslationAI_coords.fasta\
 /projects/splitorfs/work/LLMs/TranslationAI/Input/RI_transcripts_110_for_TranslationAI.fasta

head "/projects/splitorfs/work/LLMs/TranslationAI/Input/RI_transcripts_110_for_TranslationAI.fasta"

echo "Running TIS preds with TranslationAI"

# translationai -I /projects/splitorfs/work/LLMs/TranslationAI/Input/RI_transcripts_110_for_TranslationAI.fasta -t 0.0000001,0.0000001

echo "Done running TIS predictions"

mv /projects/splitorfs/work/LLMs/TranslationAI/Input/*_0.0000001.txt /projects/splitorfs/work/LLMs/TranslationAI/Output



 python /home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/analyze_TIS/analyze_TIS.py\
 /projects/splitorfs/work/LLMs/TranslationAI/Output/RI_transcripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt\
 /projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_RI.txt\
 /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample/RI_genome\
 /projects/splitorfs/work/Riboseq/Output/RiboTISH_RI_custom\
 /projects/splitorfs/work/LLMs/TranslationAI/Input/RI_transcripts_cDNA_coordinates.txt\
 /projects/splitorfs/work/LLMs/TranslationAI/Input/k4neo_val_transcripts/RI_transcripts_found_with_k4neo.txt


