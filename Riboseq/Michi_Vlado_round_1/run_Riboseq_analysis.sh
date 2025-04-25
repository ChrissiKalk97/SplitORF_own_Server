#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq

bash preprocessing_cutadapt_steps.sh
python preprocessing/cutadapt_output_parsing.py