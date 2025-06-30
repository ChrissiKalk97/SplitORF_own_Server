# BiocManager::install("Rsubread")
library(Rsubread)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)


bam_dir <- "/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"
bam_dir_2 <- "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10/NMD_genome"
ensembl_gtf <- "/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
out_dir <- "/projects/splitorfs/work/LLMs/TranslationAI/Output/analyze_TIS/genes_to_keep"



bam_files <- list.files(path = bam_dir, pattern = "filtered\\.bam$", full.names = TRUE)
bam_files_2 <- list.files(path = bam_dir_2, pattern 
= "sorted\\.bam$", full.names = TRUE)
bam_files_combined <- c(bam_files, bam_files_2)

#################################################################################
# ------------------ HELPER FUNCTIONS                        ------------------ #
#################################################################################

# TPM need function
calculate_tpm_from_feature_counts <- function(gene_counts_sample, gene_lengths){
    mapped_reads_normalized_by_length <- gene_counts_sample/gene_lengths
    tpm_values <- (mapped_reads_normalized_by_length/sum(mapped_reads_normalized_by_length)) * 1e6
    return (tpm_values)
}


gene_counts <- featureCounts(files=bam_files_combined,
annot.ext=ensembl_gtf,
isGTFAnnotationFile=TRUE,
GTF.featureType="exon", # default 
GTF.attrType="gene_id", # default
countMultiMappingReads=TRUE,
allowMultiOverlap=TRUE,
fracOverlap=0.3,
nthreads=16)


# make a df
counts_df <- as.data.frame(gene_counts$counts)
gene_lengths <- gene_counts$annotation$Length

samples <- colnames(counts_df)
for (sample_name in samples){
    gene_counts_sample <- counts_df[,sample_name]
    tpm_col_name <- paste(sample_name, "tpm", sep = "_")
    counts_df[tpm_col_name] = calculate_tpm_from_feature_counts(gene_counts_sample, gene_lengths)
}


tpm_cols <- colnames(counts_df)[str_detect(colnames(counts_df), "tpm")]
tpm_counts <- counts_df[,tpm_cols]
tpm_counts$number_above_05 <- rowSums(tpm_counts > 0.5)
tpm_filtered <- tpm_counts[tpm_counts$number_above_05 > 2,]
genes_to_keep <- rownames(tpm_filtered)



# Open PNG device first
png(file.path(out_dir, "tpm_count_genes_below_1.png"), width = 1800, height = 2400, res = 150)

# Now set layout and text formatting
par(mfrow = c(6, 3),         # 6 rows, 3 columns
    mar = c(4, 4, 2, 1),     # bottom, left, top, right margins
    cex.main = 1,            # title size
    cex.axis = 0.8,          # axis label size
    cex.lab = 0.9)           # axis title size

# Plot histograms
for (col in tpm_cols) {
  hist(tpm_counts[[col]][tpm_counts[[col]] < 1], 
       main = col, 
       xlab = "TPM", 
       breaks = seq(0, 1, by = 0.01),  # bins of width 0.01 in 0-1 range
       col = "skyblue", 
       ylim = c(0, 5000),
       xlim = c(0, 1))
}

# Close the device to write the file
dev.off()

# Open PNG device first
png(file.path(out_dir, "tpm_count_genes_greater_1.png"), width = 1800, height = 2400, res = 150)

# Now set layout and text formatting
par(mfrow = c(6, 3),         # 6 rows, 3 columns
    mar = c(4, 4, 2, 1),     # bottom, left, top, right margins
    cex.main = 1,            # title size
    cex.axis = 0.8,          # axis label size
    cex.lab = 0.9)           # axis title size

# Plot histograms
for (col in tpm_cols) {
  hist(tpm_counts[[col]][(tpm_counts[[col]] > 1) & (tpm_counts[[col]] < 1000)], 
       main = col, 
       xlab = "TPM", 
       col = "skyblue",
       xlim = c(1, 1000))
}

# Close the device to write the file
dev.off()


# Open PNG device first
png(file.path(out_dir, "tpm_count_genes.png"), width = 1800, height = 2400, res = 150)

# Now set layout and text formatting
par(mfrow = c(6, 3),         # 6 rows, 3 columns
    mar = c(4, 4, 2, 1),     # bottom, left, top, right margins
    cex.main = 1,            # title size
    cex.axis = 0.8,          # axis label size
    cex.lab = 0.9)           # axis title size

# Plot histograms
for (col in tpm_cols) {
  hist(tpm_counts[[col]][tpm_counts[[col]] < 200],
       main = col, 
       xlab = "TPM", 
       col = "skyblue", 
       breaks = 100,
       ylim = c(0, 500))
}

# Close the device to write the file
dev.off()




write(genes_to_keep, file.path(out_dir, "genes_to_keep.txt"))
