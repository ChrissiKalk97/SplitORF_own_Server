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