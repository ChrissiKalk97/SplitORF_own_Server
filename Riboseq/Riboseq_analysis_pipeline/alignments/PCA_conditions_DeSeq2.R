library(DESeq2)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]

idx_files <- c()
if (!(identical(
  list.files(in_dir,
    pattern = paste0("*", "dedup_idxstats.out")
  ),
  character(0)
))) {
  idx_files <- c(
    idx_files,
    paste0(in_dir, "/", list.files(in_dir,
      pattern = paste0("*", "dedup_idxstats.out")
    ))
  )
}

merged_counts_df <- NULL
for (file in idx_files) {
  sample <- basename(file)
  sample <- substr(sample, 20, 43)


  transcript_counts <- read.table(file,
    header = FALSE,
    sep = "\t",
    col.names = c("tID", "length", "counts", "unmapped")
  )

  transcript_counts <- transcript_counts[, c("tID", "counts")]
  if (!is.null(merged_counts_df)) {
    merged_counts_df[[sample]] <- transcript_counts$counts
  } else {
    merged_counts_df <- transcript_counts
    colnames(merged_counts_df)[2] <- sample
  }
}

rownames(merged_counts_df) <- merged_counts_df$tID
merged_counts_df <- merged_counts_df[
  ,
  !(names(merged_counts_df) %in% c("counts", "tID"))
]

merged_counts_matrix <- data.matrix(merged_counts_df)
merged_counts_matrix <- merged_counts_matrix[
  rowSums(merged_counts_matrix) > 0,
]

col_data <- data.frame(
  row.names = colnames(merged_counts_df),
  condition = factor(c("norm", "norm", "norm", "hypo", "hypo", "hypo")),
  batch = factor(c(2, 3, 4, 2, 3, 4))
)

dds <- DESeqDataSetFromMatrix(merged_counts_matrix,
  colData = col_data,
  design = ~ batch + condition
)

# DESeq() running before the PCA does not make a difference
dds <- DESeq(dds)

# this is called automatically when doing vst or rlog
vid <- vst(dds, blind = FALSE)

png(
  file.path(
    in_dir,
    "PCA_transcriptome_counts_greater_0_vid.png"
  ),
  width = 800, height = 800
)

print(plotPCA(vid, intgroup = c("condition")))
dev.off()

summary(sizeFactors(dds))

which(sizeFactors(dds) > 10)
# NULL



# RLOG TRANSFORMATION
rid <- rlog(dds, blind = FALSE)
print("correlation rid")
cor(assay(rid))

png(file.path(in_dir, "PCA_transcriptome_counts_greater_0_rlog.png"),
  width = 800, height = 800
)
print(plotPCA(rid, intgroup = c("condition")))
dev.off()


###############################################################################
# PLOT SMAPLE NAMES ON TOP
###############################################################################

pca_data <- plotPCA(rid, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_data$name <- sub("\\.+$", "", pca_data$name)

png(
  file.path(
    in_dir,
    "PCA_transcriptome_counts_greater_0_rlog_sample_names.png"
  ),
  width = 800, height = 800
)
ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1.2, size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  scale_x_continuous(limits = c(min(pca_data$PC1) - 3, max(pca_data$PC1) + 3)) +
  theme_minimal()
dev.off()




###############################################################################
# MA PLOT
###############################################################################

dds <- DESeq(dds)
res <- results(dds)
png(file.path(in_dir, "MA_plot_transcriptome.png"), width = 800, height = 800)
print(plotMA(res), ylim = c(-2, 2))
dev.off()

summary(dds)

resultsNames(dds)

res_lfc <- lfcShrink(dds, coef = 2, type = "apeglm")
png(file.path(in_dir, "MA_res_LFC_plot_transcriptome.png"),
  width = 800, height = 800
)
print(plotMA(res_lfc, ylim = c(-2, 2)))
dev.off()



norm_counts <- counts(dds, normalized = TRUE)
condition <- colData(dds)$condition
control_samples <- colnames(dds)[condition == "norm"]
print(control_samples)
treated_samples <- colnames(dds)[condition == "hypo"]
print(treated_samples)
avg_control <- rowMeans(norm_counts[, control_samples])
avg_treated <- rowMeans(norm_counts[, treated_samples])

png(file.path(in_dir, "scatter_DESeq2_all_conditions.png"),
  width = 800, height = 800
)
plot(
  x = log2(avg_control + 1),
  y = log2(avg_treated + 1),
  xlab = "log2(Average norm expression + 1)",
  ylab = "log2(Average hypo expression + 1)",
  main = "Gene Expression: hypo vs norm",
  pch = 16, col = rgb(0, 0, 1, 0.3)
)
abline(0, 1, col = "red") # line of equality
dev.off()


################################################################################
################################################################################
################################################################################


################################################################################
# FILTER FOR mRNA mapping reads                                                #
################################################################################
mane_mrnas <- read.table(
  "/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai", # nolint: line_length_linter.
  header = FALSE,
  sep = "\t",
  col.names = c("tID", "length", "some_nr", "some_other_nr", "yet_another_nr")
) # nolint: line_length_linter.

merged_counts_df_mrna <- merged_counts_df[
  rownames(merged_counts_df) %in% mane_mrnas$tID,
]

merged_counts_matrix_mrna <- data.matrix(merged_counts_df_mrna)
merged_counts_matrix_mrna <- merged_counts_matrix_mrna[
  rowSums(merged_counts_matrix_mrna) > 0,
]

col_data <- data.frame(
  row.names = colnames(merged_counts_df_mrna),
  condition = factor(c("norm", "norm", "norm", "hypo", "hypo", "hypo")),
  batch = factor(c(2, 3, 4, 2, 3, 4))
)

dds_mrna <- DESeqDataSetFromMatrix(merged_counts_matrix_mrna,
  colData = col_data, design = ~ batch + condition
)

# DESeq() running before the PCA does not make a difference
dds_mrna <- DESeq(dds_mrna)
# this is called automatically when doing vst or rlog
vid_mrna <- vst(dds_mrna, blind = FALSE)

png(file.path(in_dir, "PCA_transcriptome_counts_greater_0_vid_mRNA.png"),
  width = 800, height = 800
)
print(plotPCA(vid_mrna, intgroup = c("condition")))
dev.off()


# RLOG TRANSFORMATION
rid_mrna <- rlog(dds_mrna, blind = FALSE)
print("correlation rid")
cor(assay(rid_mrna))

png(file.path(in_dir, "PCA_transcriptome_counts_greater_0_rlog_mRNA.png"),
  width = 800, height = 800
)
print(plotPCA(rid_mrna, intgroup = c("condition")))
dev.off()


###############################################################################
# PLOT SMAPLE NAMES ON TOP
###############################################################################

pca_data_mrna <- plotPCA(rid_mrna, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data_mrna, "percentVar"))
pca_data_mrna$name <- sub("\\.+$", "", pca_data_mrna$name)

png(
  file.path(
    in_dir,
    "PCA_transcriptome_counts_greater_0_rlog_sample_names_mRNA.png"
  ),
  width = 800, height = 800
)
ggplot(pca_data_mrna, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1.2, size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  scale_x_continuous(limits = c(
    min(pca_data_mrna$PC1) - 3,
    max(pca_data_mrna$PC1) + 3
  )) +
  theme_minimal()
dev.off()




###############################################################################
# MA PLOT
###############################################################################

dds_mrna <- DESeq(dds_mrna)
res_mrna <- results(dds_mrna, contrast = c("condition", "hypo", "norm"))
png(file.path(in_dir, "MA_plot_transcriptome_mRNA.png"),
  width = 800, height = 800
)
print(plotMA(res_mrna), ylim = c(-2, 2))
dev.off()

summary(dds_mrna)

resultsNames(dds_mrna)

res_lfc_mrna <- lfcShrink(dds_mrna, coef = 2, type = "apeglm")
png(file.path(in_dir, "MA_res_LFC_plot_transcriptome_mRNA.png"),
  width = 800, height = 800
)
print(plotMA(res_lfc_mrna, ylim = c(-2, 2)))
dev.off()



norm_counts <- counts(dds_mrna, normalized = TRUE)
condition <- colData(dds_mrna)$condition
control_samples_mrna <- colnames(dds_mrna)[condition == "norm"]
print(control_samples_mrna)
treated_samples_mrna <- colnames(dds_mrna)[condition == "hypo"]
print(treated_samples_mrna)
avg_control_mrna <- rowMeans(norm_counts[, control_samples_mrna])
avg_treated_mrna <- rowMeans(norm_counts[, treated_samples_mrna])

png(file.path(in_dir, "scatter_DESeq2_all_conditions_mRNA.png"),
  width = 800, height = 800
)
plot(
  x = log2(avg_control_mrna + 1),
  y = log2(avg_treated_mrna + 1),
  xlab = "log2(Average norm expression + 1)",
  ylab = "log2(Average hypo expression + 1)",
  main = "Gene Expression: hypo vs norm",
  pch = 16, col = rgb(0, 0, 1, 0.3)
)
abline(0, 1, col = "red") # line of equality
dev.off()


###############################################################################
# Differentially expressed genes
###############################################################################
res_mrna <- res_mrna[res_mrna$baseMean > 0, ]
res_mrna <- na.omit(res_mrna)
res_mrna_pval <- res_mrna[res_mrna$pvalue < 0.05, ]

# write output
writeLines(
  rownames(res_mrna_pval[
    (res_mrna_pval$padj < 0.05) & (res_mrna_pval$log2FoldChange > 1),
  ]),
  file.path(in_dir, "DEGs", "hypo_vs_norm_Riboseq_DEGs_1.txt")
)


writeLines(
  rownames(res_mrna_pval[
    (res_mrna_pval$padj < 0.05) & (res_mrna_pval$log2FoldChange < -1),
  ]),
  file.path(in_dir, "DEGs", "hypo_vs_norm_downreg_Riboseq_DEGs_1.txt")
)


writeLines(
  rownames(res_mrna_pval[
    (res_mrna_pval$padj < 0.05) & ((res_mrna_pval$log2FoldChange < -1) | (res_mrna_pval$log2FoldChange > 1)),
  ]),
  file.path(in_dir, "DEGs", "hypo_vs_norm_up_down_Riboseq_DEGs_1.txt")
)

