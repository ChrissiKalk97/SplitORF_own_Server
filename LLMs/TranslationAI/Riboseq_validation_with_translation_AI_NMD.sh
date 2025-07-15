#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate translationAI

# cd /projects/splitorfs/work/LLMs/TranslationAI/Output

# translationai -I /projects/splitorfs/work/LLMs/TranslationAI/Input/NMD_trnascripts_110_for_TranslationAI.fasta -t 0.0000001,0.0000001

# python /home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/analyze_TIS/analyze_TIS.py\
#  /projects/splitorfs/work/LLMs/TranslationAI/Output/NMD_trnascripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt\
#  /projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt\
#  /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample/NMD_genome\
#  /projects/splitorfs/work/Riboseq/Output/RiboTISH_NMD_custom\
#  /projects/splitorfs/work/LLMs/TranslationAI/Input/NMD_transcripts_cDNA_coordinates.txt\
#  /projects/splitorfs/work/LLMs/TranslationAI/Input/k4neo_val_transcripts/NMD_transcripts_found_with_k4neo.txt


# Rscript /home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/analyze_TIS/filter_genes_by_tpm.R

python /home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/analyze_TIS/analyze_TIS_pipeline.py \
  /projects/splitorfs/work/LLMs/TranslationAI/Output/NMD_trnascripts_110_for_TranslationAI.fa_predTIS_0.0000001.txt \
  /projects/splitorfs/work/LLMs/TIS_transformer/Input/SO_pipeline_results/UniqueProteinORFPairs_NMD.txt \
  /projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/NMD_genome \
  /projects/splitorfs/work/Riboseq/Output/RiboTISH_NMD_custom \
  --genes_to_keep_file /projects/splitorfs/work/LLMs/TranslationAI/Output/analyze_TIS/genes_to_keep/genes_to_keep.txt \
  --k4neo_path /projects/splitorfs/work/LLMs/TranslationAI/Input/NMD_transcripts_cDNA_coordinates.txt \
  --Ensembl_canonical_path /projects/splitorfs/work/LLMs/TranslationAI/Input/k4neo_val_transcripts/NMD_transcripts_found_with_k4neo.txt


