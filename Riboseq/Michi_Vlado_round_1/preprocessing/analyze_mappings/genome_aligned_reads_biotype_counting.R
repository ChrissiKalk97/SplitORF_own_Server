# BiocManager::install("Rsubread")
library(Rsubread)
library(tidyr)
library(dplyr)


try <- featureCounts(files="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated/uf_muellermcnicoll_2025_04_01_huvec_dnor_2_dedup_filtered.bam",
annot.ext="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf",
# annot.ext="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf",
isGTFAnnotationFile=TRUE,
GTF.featureType="gene",
GTF.attrType="gene_biotype",
countMultiMappingReads=TRUE,
fraction=TRUE,
allowMultiOverlap=TRUE,
fracOverlap=0.3,
nthreads=16)

try$stat


counts <- try$counts
anno <- try$annotation

# Combine biotype and counts
# biotype_counts <- data.frame(
#   biotype = anno$GeneID,  # This will be the biotype if GTF.attrType="gene_biotype"
#   counts = counts[, 1]    # Assumes one sample/BAM
# )

# # Sum counts by biotype
# biotype_summary <- aggregate(counts ~ biotype, data = biotype_counts, FUN = sum)

# # View result
# biotype_summary[order(-biotype_summary$counts), ]