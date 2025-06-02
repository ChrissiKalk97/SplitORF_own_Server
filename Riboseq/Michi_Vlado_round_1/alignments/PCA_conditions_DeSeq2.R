library(DESeq2)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- args[1]

idx_files = c()
if(!(identical(list.files(in_dir, pattern=paste0("*","dedup_idxstats.out")), character(0)))){
      idx_files <- c(idx_files, paste0(in_dir,"/",list.files(in_dir, pattern=paste0("*","dedup_idxstats.out"))))
    }

merged_counts_df <- NULL
for (file in idx_files){
    sample <- basename(file)
    sample <- substr(sample, 20, 43)


    transcript_counts <- read.table(file,
    header = FALSE,
    sep = "\t",
    col.names = c('tID', 'length', 'counts', 'unmapped'))

    transcript_counts <- transcript_counts[,c('tID', 'counts')]
    if (!is.null(merged_counts_df)) {
        merged_counts_df[[sample]] <- transcript_counts$counts
    }
    else {
        merged_counts_df <- transcript_counts
        colnames(merged_counts_df)[2] <- sample
        

        
    }

}

rownames(merged_counts_df) <- merged_counts_df$tID
merged_counts_df <- merged_counts_df[ , !(names(merged_counts_df) %in% c('counts', 'tID'))]

merged_counts_matrix <- data.matrix(merged_counts_df)
merged_counts_matrix <- merged_counts_matrix[rowSums(merged_counts_matrix) > 0,]

col_data <- data.frame(row.names = colnames(merged_counts_df),
 condition = c('norm', 'norm', 'norm', 'hypo', 'hypo', 'hypo'))

dds <- DESeqDataSetFromMatrix(merged_counts_matrix, colData = col_data, design = ~condition)

# DESeq() running before the PCA does not make a difference
dds <- DESeq(dds)
# dds <- estimateSizeFactors(dds)
# this is called automatically when doing vst or rlog
vid <- vst(dds, blind = FALSE)

png(file.path(in_dir, 'PCA_transcriptome_counts_greater_0_vid.png'), width = 800, height = 800)
print(plotPCA(vid, intgroup = c("condition")))
dev.off()

summary(sizeFactors(dds))

which(sizeFactors(dds) > 10)
# NULL



# RLOG TRANSFORMATION
rid <- rlog(dds, blind = FALSE)
print('correlation rid')
cor(assay(rid))

png(file.path(in_dir, 'PCA_transcriptome_counts_greater_0_rlog.png'), width = 800, height = 800)
print(plotPCA(rid, intgroup = c("condition")))
dev.off()


###############################################################################
#PLOT SMAPLE NAMES ON TOP
###############################################################################

pcaData <- plotPCA(rid, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$name <- sub("\\.+$", "", pcaData$name)

png(file.path(in_dir, 'PCA_transcriptome_counts_greater_0_rlog_sample_names.png'), width = 800, height = 800)
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1.2, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_x_continuous(limits = c(min(pcaData$PC1) - 3, max(pcaData$PC1) + 3)) +
  theme_minimal()
dev.off()




###############################################################################
# MA PLOT
###############################################################################

dds <- DESeq(dds)
res <- results(dds)
png(file.path(in_dir, 'MA_plot_transcriptome.png'), width = 800, height = 800)
print(plotMA(res), ylim=c(-2,2))
dev.off()

summary(dds)

resultsNames(dds)

resLFC <- lfcShrink(dds, coef = 2, type="apeglm")
png(file.path(in_dir, 'MA_res_LFC_plot_transcriptome.png'), width = 800, height = 800)
print(plotMA(resLFC, ylim=c(-2,2)))
dev.off()



norm_counts <- counts(dds, normalized = TRUE)
condition <- colData(dds)$condition
control_samples <- colnames(dds)[condition == "norm"]
print(control_samples)
treated_samples <- colnames(dds)[condition == "hypo"]
print(treated_samples)
avg_control <- rowMeans(norm_counts[, control_samples])
avg_treated <- rowMeans(norm_counts[, treated_samples])

png(file.path(in_dir, 'scatter_DESeq2_all_conditions.png'), width = 800, height = 800)
plot(
  x = log2(avg_control + 1), 
  y = log2(avg_treated + 1),
  xlab = "log2(Average norm expression + 1)",
  ylab = "log2(Average hypo expression + 1)",
  main = "Gene Expression: hypo vs norm",
  pch = 16, col = rgb(0, 0, 1, 0.3)
)
abline(0, 1, col = "red")  # line of equality
dev.off()


################################################################################
################################################################################
################################################################################


################################################################################
# FILTER FOR mRNA mapping reads                                                #
################################################################################
MANE_mRNAs <- read.table('/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai',
    header = FALSE,
    sep = "\t",
    col.names = c('tID', 'length', 'some_nr', 'some_other_nr', 'yet_another_nr'))

merged_counts_df_mRNA <- merged_counts_df[rownames(merged_counts_df) %in% MANE_mRNAs$tID, ]

merged_counts_matrix_mRNA <- data.matrix(merged_counts_df_mRNA)
merged_counts_matrix_mRNA <- merged_counts_matrix_mRNA[rowSums(merged_counts_matrix_mRNA) > 0,]

col_data <- data.frame(row.names = colnames(merged_counts_df_mRNA),
 condition = c('norm', 'norm', 'norm', 'hypo', 'hypo', 'hypo'))

dds_mRNA <- DESeqDataSetFromMatrix(merged_counts_matrix_mRNA, colData = col_data, design = ~condition)

# DESeq() running before the PCA does not make a difference
dds_mRNA <- DESeq(dds_mRNA)
# dds_mRNA <- estimateSizeFactors(dds_mRNA)
# this is called automatically when doing vst or rlog
vid_mRNA <- vst(dds_mRNA, blind = FALSE)

png(file.path(in_dir, 'PCA_transcriptome_counts_greater_0_vid_mRNA.png'), width = 800, height = 800)
print(plotPCA(vid_mRNA, intgroup = c("condition")))
dev.off()

summary(sizeFactors(dds_mRNA))

which(sizeFactors(dds_mRNA) > 10)
# NULL



# RLOG TRANSFORMATION
rid_mRNA <- rlog(dds_mRNA, blind = FALSE)
print('correlation rid')
cor(assay(rid_mRNA))

png(file.path(in_dir, 'PCA_transcriptome_counts_greater_0_rlog_mRNA.png'), width = 800, height = 800)
print(plotPCA(rid_mRNA, intgroup = c("condition")))
dev.off()


###############################################################################
# PLOT SMAPLE NAMES ON TOP
###############################################################################

pcaData_mRNA <- plotPCA(rid_mRNA, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData_mRNA, "percentVar"))
pcaData_mRNA$name <- sub("\\.+$", "", pcaData_mRNA$name)

png(file.path(in_dir, 'PCA_transcriptome_counts_greater_0_rlog_sample_names_mRNA.png'), width = 800, height = 800)
ggplot(pcaData_mRNA, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1.2, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_x_continuous(limits = c(min(pcaData_mRNA$PC1) - 3, max(pcaData_mRNA$PC1) + 3)) +
  theme_minimal()
dev.off()




###############################################################################
# MA PLOT
###############################################################################

dds_mRNA <- DESeq(dds_mRNA)
res_mRNA <- results(dds_mRNA)
png(file.path(in_dir, 'MA_plot_transcriptome_mRNA.png'), width = 800, height = 800)
print(plotMA(res_mRNA), ylim=c(-2,2))
dev.off()

summary(dds_mRNA)

resultsNames(dds_mRNA)

resLFC_mRNA <- lfcShrink(dds_mRNA, coef = 2, type="apeglm")
png(file.path(in_dir, 'MA_res_LFC_plot_transcriptome_mRNA.png'), width = 800, height = 800)
print(plotMA(resLFC_mRNA, ylim=c(-2,2)))
dev.off()



norm_counts <- counts(dds_mRNA, normalized = TRUE)
condition <- colData(dds_mRNA)$condition
control_samples_mRNA <- colnames(dds_mRNA)[condition == "norm"]
print(control_samples_mRNA)
treated_samples_mRNA <- colnames(dds_mRNA)[condition == "hypo"]
print(treated_samples_mRNA)
avg_control_mRNA <- rowMeans(norm_counts[, control_samples_mRNA])
avg_treated_mRNA <- rowMeans(norm_counts[, treated_samples_mRNA])

png(file.path(in_dir, 'scatter_DESeq2_all_conditions_mRNA.png'), width = 800, height = 800)
plot(
  x = log2(avg_control_mRNA + 1), 
  y = log2(avg_treated_mRNA + 1),
  xlab = "log2(Average norm expression + 1)",
  ylab = "log2(Average hypo expression + 1)",
  main = "Gene Expression: hypo vs norm",
  pch = 16, col = rgb(0, 0, 1, 0.3)
)
abline(0, 1, col = "red")  # line of equality
dev.off()
