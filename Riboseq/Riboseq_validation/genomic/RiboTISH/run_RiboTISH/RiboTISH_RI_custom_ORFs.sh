#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq
#This script calls the aligning scripts to align the Ribo-Seq reads to the transcriptome
# and calculates the overlap with the previously determined unique regions (main pipeline) of RI and NMD
#The following are the main directories for the pipeline outputs and for the alignment outputs
genomeFasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
genomeGTF="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
outputRiboTISH="/projects/splitorfs/work/Riboseq/Output/RiboTISH_RI_custom"
outputRiboTISH_old="/projects/splitorfs/work/Riboseq/Output/RiboTISH"


if [ ! -d $outputRiboTISH ];then
	mkdir $outputRiboTISH
fi


# Input Riboseq fastq files
input_data="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample/RI_genome/*_sorted.bam"



# ribotish quality -b $input_data -g $genomeGTF -f $outputRiboTISH/ERR3367797_qual.pdf -r $outputRiboTISH/ERR3367797_para.py -o $outputRiboTISH/ERR3367797_qual.txt

# ribotish predict -b $input_data -g $genomeGTF -f $genomeFasta --longest -o $outputRiboTISH/ERR3367797_pred_longest.txt



# ribotish predict -b $input_data -g $genomeGTF -f $genomeFasta --framebest -o $outputRiboTISH/ERR3367797_pred_framebest.txt

# -i INPUT

# Only test input candidate ORFs, format: transcript ID, start, stop

# This is very handy for Split-ORFs!


for i in $input_data; do
	sample_name=$(basename "$i" _fastp.fastq)
	# ribotish quality\
	# -b $i\
	# -g $genomeGTF\
	# -f $outputRiboTISH/${sample_name}_qual.pdf\
	# -r $outputRiboTISH/${sample_name}_para.py\
	# -o $outputRiboTISH/${sample_name}_qual.txt

	ribotish predict\
	-b $i\
	-g $genomeGTF\
	-f $genomeFasta\
	--framebest \
	--inframecount\
	--ribopara $outputRiboTISH_old/${sample_name}_para.py\
    -i /projects/splitorfs/work/Riboseq/data/RiboTISH/Unique_DNA_Regions_for_riboseq_RI.bed_RiboTISH_modified.bed\
    --transprofile $outputRiboTISH/${sample_name}_transprofile.py\
	-o $outputRiboTISH/${sample_name}_pred_framebest.txt


	echo "Predicted translation for sample ${sample_name}"
done

python analyze_predicted_URs_riboTISH.py /projects/splitorfs/work/Riboseq/Output/RiboTISH_RI_custom

python compare_URs_riboTISH_with_emp_background.py /projects/splitorfs/work/Riboseq/Output/RiboTISH_RI_custom
