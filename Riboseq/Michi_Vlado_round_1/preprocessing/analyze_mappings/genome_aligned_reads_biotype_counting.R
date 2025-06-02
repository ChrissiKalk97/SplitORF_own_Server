# BiocManager::install("Rsubread")
library(Rsubread)
library(tidyr)
library(dplyr)
library(ggplot2)


bam_dir <- "/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"
bam_files <- list.files(path = bam_dir, pattern = "filtered\\.bam$", full.names = TRUE)

try <- featureCounts(files=bam_files,
annot.ext="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf",
isGTFAnnotationFile=TRUE,
GTF.featureType="gene",
GTF.attrType="gene_biotype",
countMultiMappingReads=TRUE,
allowMultiOverlap=TRUE,
fracOverlap=0.3,
nthreads=16)


counts <- try$counts

# make a df
counts_df <- as.data.frame(counts)
# only keep relevant counts
# counts_df <- counts_df[rowSums(counts_df) > 20,]
# biotype in index, but need as column
counts_df$biotype <- rownames(counts_df)

counts_df_longer <- pivot_longer(counts_df, cols = colnames(counts_df)[1:6], names_to = "sample")

counts_df_longer$biotypes_combined <- "other_ncRNA"
counts_df_longer[counts_df_longer$biotype == 'protein_coding', 'biotypes_combined'] <- 'protein_coding'
counts_df_longer[counts_df_longer$biotype == 'lncRNA', 'biotypes_combined'] <- 'lncRNA'
counts_df_longer[counts_df_longer$biotype == 'rRNA', 'biotypes_combined'] <- 'rRNA'
counts_df_longer[counts_df_longer$biotype == 'rRNA_pseudogene', 'biotypes_combined'] <- 'rRNA_pseudogene'
counts_df_longer[counts_df_longer$biotype == 'Mt_tRNA', 'biotypes_combined'] <- 'Mt_tRNA'
counts_df_longer[counts_df_longer$biotype == 'Mt_rRNA', 'biotypes_combined'] <- 'Mt_rRNA'



# geom_col has stat_identity as unchangeable default
# stat identity means we provide the values as y, not counting!
barplot_stacked <- ggplot(counts_df_longer, aes(fill=biotypes_combined, y=value, x=sample)) + 
    geom_col(position="fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(bam_dir, "barplot_mapping_stats_biotypes_combined.svg"), plot = barplot_stacked, width = 12, height = 12, units = "in", dpi = 600)

barplot_stacked <- ggplot(counts_df_longer, aes(fill=biotypes_combined, y=value, x=sample)) + 
  geom_col(position="fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(bam_dir, "barplot_mapping_stats_biotypes_combined.png"), plot = barplot_stacked, width = 12, height = 12, units = "in", dpi = 600)

barplot_stacked <- ggplot(counts_df_longer, aes(fill=biotypes_combined, y=value, x=sample)) + 
  geom_col(position="fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(bam_dir, "barplot_mapping_stats_biotypes_combined.pdf"), plot = barplot_stacked, width = 12, height = 12, units = "in", dpi = 600)

barplot_stacked <- ggplot(counts_df_longer, aes(fill=biotype, y=value, x=sample)) + 
    geom_col(position="fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(bam_dir, "barplot_mapping_stats.svg"), plot = barplot_stacked, width = 12, height = 12, units = "in", dpi = 600)