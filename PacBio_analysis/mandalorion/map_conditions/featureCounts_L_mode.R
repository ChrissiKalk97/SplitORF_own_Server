library(Rsubread)
library(tidyr)
library(dplyr)
library(ggplot2)




args <- commandArgs(trailingOnly = TRUE)
bam_dir <- args[1]
gtf <- args[2]

bam_files <- list.files(path = bam_dir, pattern = "\\.bam$", full.names = TRUE)

feature_counts_long <- featureCounts(files=bam_files,
annot.ext=gtf,
isGTFAnnotationFile=TRUE,
GTF.featureType="transcript",
countMultiMappingReads=TRUE,
allowMultiOverlap=TRUE,
fracOverlap=0.8,
isLongRead=TRUE,
nthreads=16)


print(feature_counts_long)
