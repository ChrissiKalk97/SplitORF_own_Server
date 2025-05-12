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
 condition = c('norm', 'norm', 'norm', 'hypo', 'hypo', 'hypo'))

dds <- DESeqDataSetFromMatrix(merged_counts_matrix, colData = col_data, design = ~condition)
# dds <- estimateSizeFactors(dds)
# this is called automatically when doing vst or rlog
vid <- vst(dds, blind = FALSE)

png(file.path(in_dir, 'PCA_transcriptome_counts_greater_0_vid.png'), width = 800, height = 800)
print(plotPCA(vid, intgroup = c("condition")))
dev.off()

summary(sizeFactors(dds))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.2796  0.3392  0.4362 16.0354  0.6618 94.0167 
which(sizeFactors(dds) > 10)
# 2025_04_01_huvec_dnor_2. 
#                        1 



# RLOG TRANSFORMATION
rid <- rlog(dds, blind = FALSE)

png(file.path(in_dir, 'PCA_transcriptome_counts_greater_0_rlog.png'), width = 800, height = 800)
print(plotPCA(rid, intgroup = c("condition")))
dev.off()
