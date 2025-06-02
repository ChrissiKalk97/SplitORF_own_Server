library(DESeq2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- args[1]

idx_files = c()
if(!(identical(list.files(in_dir, pattern=paste0("*","idxstats.out")), character(0)))){
      idx_files <- c(idx_files, paste0(in_dir,"/",list.files(in_dir, pattern=paste0("*","idxstats.out"))))
    }

merged_counts_df <- NULL
for (file in idx_files){
    sample <- basename(file)
    sample <- substr(sample, 20, 43)
    sample <- str_split(sample, '\\.')[[1]][1]

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
 condition = c('SeRP_CHX_A2', 
              'SeRP_CHX_A2', 
              'SeRP_CHX_A2',
              'SeRP_CHX_A2',
              'SeRP_Puro_A2', 
              'SeRP_Puro_A2', 
              'SeRP_Puro_A2',
              'SeRP_CHX_B1', 
              'SeRP_CHX_B1', 
              'SeRP_CHX_B1', 
              'SeRP_CHX_B1', 
              'SeRP_Puro_B1', 
              'SeRP_Puro_B1', 
              'SeRP_Puro_B1', 
              'In_CHX',
              'In_CHX', 
              'In_CHX',
              'In_CHX', 
              'In_Puro',
              'In_Puro', 
              'In_Puro', 
              'SeRP_CHX_M', 
              'SeRP_CHX_M',
              'SeRP_CHX_M', 
              'SeRP_CHX_M'
 ))

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



################################################################################
################################################################################
################################################################################


################################################################################
# OBTAINING DEGs (transcripts) between SeRP and Input ##########################
################################################################################

# SERP A2 CHX
dds_mRNA_CHX_A2 <- dds_mRNA[,dds_mRNA$condition %in% c('SeRP_CHX_A2', 'In_CHX')]
dds_mRNA_CHX_A2 <- dds_mRNA_CHX_A2[rowSums(counts(dds_mRNA_CHX_A2)) >= 10, ]
res_CHX_A2 <- results(dds_mRNA_CHX_A2, contrast=c("condition", "SeRP_CHX_A2","In_CHX"))
# kept the independent filter, but this does not change the number of diff genes
res_CHX_A2 <-  res_CHX_A2[res_CHX_A2$baseMean > 0, ]
res_CHX_A2_pval <- res_CHX_A2[res_CHX_A2$pvalue < 0.05, ]
# rownames(res_CHX_A2_pval[res_CHX_A2_pval$padj < 0.05, ])
writeLines(rownames(res_CHX_A2[res_CHX_A2$padj < 0.05, ]), file.path(in_dir, "DEGs_CHX_A2_mRNA.txt"))

# SERP B1 CHX
dds_mRNA_CHX_B1 <- dds_mRNA[,dds_mRNA$condition %in% c('SeRP_CHX_B1', 'In_CHX')]
dds_mRNA_CHX_B1 <- dds_mRNA_CHX_B1[rowSums(counts(dds_mRNA_CHX_B1)) >= 10, ]
res_CHX_B1 <- results(dds_mRNA_CHX_B1, contrast=c("condition", "SeRP_CHX_B1","In_CHX"))
res_CHX_B1 <-  res_CHX_B1[res_CHX_B1$baseMean > 0, ]
res_CHX_B1_pval <- res_CHX_B1[res_CHX_B1$pvalue < 0.05, ]
# rownames(res_CHX_B1_pval[res_CHX_B1_pval$padj < 0.05, ])
writeLines(rownames(res_CHX_B1[res_CHX_B1$padj < 0.05, ]), file.path(in_dir, "DEGs_CHX_B1_mRNA.txt"))

# SERP MOCK CHX
dds_mRNA_CHX_M <- dds_mRNA[,dds_mRNA$condition %in% c('SeRP_CHX_M', 'In_CHX')]
dds_mRNA_CHX_M <- dds_mRNA_CHX_M[rowSums(counts(dds_mRNA_CHX_M)) >= 10, ]
res_CHX_M <- results(dds_mRNA_CHX_M,contrast=c("condition", "SeRP_CHX_M","In_CHX"))
res_CHX_M <-  res_CHX_M[res_CHX_M$baseMean > 0, ]
res_CHX_M_pval <- res_CHX_M[res_CHX_M$pvalue < 0.05, ]
# rownames(res_CHX_M_pval[res_CHX_M_pval$padj < 0.05, ])
writeLines(rownames(res_CHX_M[res_CHX_M$padj < 0.05, ]), file.path(in_dir, "DEGs_CHX_M_mRNA.txt"))

# SERP A2 Puro
# this reaises error: 
# res_Puro_A2 <- results(dds_mRNA_Puro_A2, contrast=c("condition", "SeRP_Puro_A2","In_Puro"))
# Error in counts(object) %*% whichSamples : non-conformable arguments
# dds_mRNA_Puro_A2 <- dds_mRNA[,dds_mRNA$condition %in% c('SeRP_Puro_A2', 'In_Puro')]
# dds_mRNA_Puro_A2 <- dds_mRNA_Puro_A2[rowSums(counts(dds_mRNA_Puro_A2)) >= 10, ]
# dds_mRNA_Puro_A2 <- DESeq(dds_mRNA_Puro_A2)
# res_Puro_A2 <- results(dds_mRNA_Puro_A2, contrast=c("condition", "SeRP_Puro_A2","In_Puro"))
# res_Puro_A2 <-  res_Puro_A2[res_Puro_A2$baseMean > 0, ]
# res_Puro_A2_pval <- res_Puro_A2[res_Puro_A2$pvalue < 0.05, ]
# # rownames(res_Puro_A2_pval[res_Puro_A2_pval$padj < 0.05, ])
# writeLines(rownames(res_Puro_A2[res_Puro_A2$padj < 0.05, ]), file.path(in_dir, "DEGs_Puro_A2_mRNA.txt"))