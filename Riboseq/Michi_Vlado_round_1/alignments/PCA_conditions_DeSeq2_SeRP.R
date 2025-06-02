library(DESeq2)

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
 condition = c('In_Puro', 'In_Puro', 'In_Puro', 'In_CHX', 'In_CHX', 'In_CHX', 
 'IP_Puro', 'IP_Puro', 'IP_Puro', 'IP_CHX', 'IP_CHX', 'IP_CHX'))

dds <- DESeqDataSetFromMatrix(merged_counts_matrix, colData = col_data, design = ~condition)
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



# MA plot
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
control_samples <- colnames(dds)[condition == "In_Puro"]
print(control_samples)
treated_samples <- colnames(dds)[condition == "IP_Puro"]
print(treated_samples)
avg_control <- rowMeans(norm_counts[, control_samples])
avg_treated <- rowMeans(norm_counts[, treated_samples])

png(file.path(in_dir, 'scatter_DESeq2_Puro.png'), width = 800, height = 800)
plot(
  x = log2(avg_control + 1), 
  y = log2(avg_treated + 1),
  xlab = "log2(Average In_Puro expression + 1)",
  ylab = "log2(Average SeRP_Puro expression + 1)",
  main = "Gene Expression: In Puro vs SeRP Puro",
  pch = 16, col = rgb(0, 0, 1, 0.3)
)
abline(0, 1, col = "red")  # line of equality
dev.off()



control_samples <- colnames(dds)[condition == "In_CHX"]
print(control_samples)
treated_samples <- colnames(dds)[condition == "IP_CHX"]
print(treated_samples)
avg_control <- rowMeans(norm_counts[, control_samples])
avg_treated <- rowMeans(norm_counts[, treated_samples])

png(file.path(in_dir, 'scatter_DESeq2_CHX.png'), width = 800, height = 800)
plot(
  x = log2(avg_control + 1), 
  y = log2(avg_treated + 1),
  xlab = "log2(Average In_CHX expression + 1)",
  ylab = "log2(Average SeRP_CHX expression + 1)",
  main = "Gene Expression: In CHX vs SeRP Puro",
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
 condition = c('In_Puro', 'In_Puro', 'In_Puro', 'In_CHX', 'In_CHX', 'In_CHX', 
 'IP_Puro', 'IP_Puro', 'IP_Puro', 'IP_CHX', 'IP_CHX', 'IP_CHX'))

dds_mRNA <- DESeqDataSetFromMatrix(merged_counts_df_mRNA, colData = col_data, design = ~condition)
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



# MA plot
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



norm_counts_mRNA <- counts(dds_mRNA, normalized = TRUE)
condition <- colData(dds_mRNA)$condition
control_samples_mRNA <- colnames(dds_mRNA)[condition == "In_Puro"]
treated_samples_mRNA <- colnames(dds_mRNA)[condition == "IP_Puro"]
avg_control_mRNA <- rowMeans(norm_counts_mRNA[, control_samples_mRNA])
avg_treated_mRNA <- rowMeans(norm_counts_mRNA[, treated_samples_mRNA])

png(file.path(in_dir, 'scatter_DESeq2_Puro_mRNA.png'), width = 800, height = 800)
plot(
  x = log2(avg_control_mRNA + 1), 
  y = log2(avg_treated_mRNA + 1),
  xlab = "log2(Average In_Puro expression + 1)",
  ylab = "log2(Average SeRP_Puro expression + 1)",
  main = "Gene Expression: In Puro vs SeRP Puro",
  pch = 16, col = rgb(0, 0, 1, 0.3)
)
abline(0, 1, col = "red")  # line of equality
dev.off()



control_samples_mRNA <- colnames(dds_mRNA)[condition == "In_CHX"]
treated_samples_mRNA <- colnames(dds_mRNA)[condition == "IP_CHX"]
avg_control_mRNA <- rowMeans(norm_counts_mRNA[, control_samples_mRNA])
avg_treated_mRNA <- rowMeans(norm_counts_mRNA[, treated_samples_mRNA])

png(file.path(in_dir, 'scatter_DESeq2_CHX_mRNA.png'), width = 800, height = 800)
plot(
  x = log2(avg_control_mRNA + 1), 
  y = log2(avg_treated_mRNA + 1),
  xlab = "log2(Average In_CHX expression + 1)",
  ylab = "log2(Average SeRP_CHX expression + 1)",
  main = "Gene Expression: In CHX vs SeRP Puro",
  pch = 16, col = rgb(0, 0, 1, 0.3)
)
abline(0, 1, col = "red")  # line of equality
dev.off()






################################################################################
################################################################################
################################################################################


################################################################################
# OBTAINING DEGs (transcripts) between SeRP and Input ##########################
################################################################################

dds_mRNA <- DESeq(dds_mRNA)
res_CHX <- results(dds_mRNA, , contrast=c("condition", "In_CHX","IP_CHX"), independentFiltering = FALSE)
# needed to switch the independent filter off: a lot of NAs were introduced
# independent filter removes genes with too low overall mean
res_CHX <-  res_CHX[res_CHX$baseMean > 0, ]
rownames(res_CHX[res_CHX$padj < 0.05, ])

writeLines(rownames(res_CHX[res_CHX$padj < 0.05, ]), file.path(in_dir, "DEGs_CHX_mRNA.txt"))