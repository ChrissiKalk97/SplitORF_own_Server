# analyze_soft_clipping.py obtains BAM files as input and counts how
# often and how many bases were soft-clipped

# uasge: python analyze_soft_clipping.py bam_directory

import sys
import os.path
import pysam
import re


bam_dir = sys.argv[1]
bam_dir = "/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1"
bam_dir = "/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia"


# 
# "2025_04_01_huvec_dnor_2.cutadapt_umi_fastp.bowtie2_concat_transcriptome_k1_sorted_filtered.bam"

soft_clip_counter_front = 0
soft_clip_counter_end = 0
soft_clip_counter = 0
trimmed_A_at_end = 0
for filename in os.listdir(bam_dir):
    if filename.endswith("2025_04_01_huvec_dnor_2.cutadapt_umi_fastp.only_R1_Aligned.sortedByCoord.out.bam"):
        bamfile = pysam.AlignmentFile(os.path.join(bam_dir, filename), "rb")
        for aln in bamfile:
            if 'S' in set(aln.cigarstring):
                groups_front = re.search(r'^[0-9]+S', aln.cigarstring)
                groups_end = re.search(r'[0-9]+S$', aln.cigarstring)
                if groups_front is not None:
                    soft_clip_counter_front += 1
                if groups_end is not None:
                    soft_clip_counter_end += 1
                    if aln.query_sequence.endswith('A'):
                        trimmed_A_at_end += 1
                soft_clip_counter += 1
