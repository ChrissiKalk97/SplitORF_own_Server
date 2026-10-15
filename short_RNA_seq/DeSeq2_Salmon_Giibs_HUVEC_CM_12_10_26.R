# BiocManager::install("edgeR")
library(edgeR)
# BiocManager::install("tximport")
library(tximport)
library(stringr)

# if (!requireNamespace("tidyverse", quietly = TRUE)) {
#     install.packages("tidyverse", repos = "https://cloud.r-project.org/")
# }
library(tidyverse)
# BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
# BiocManager::install("txdbmaker")
library(txdbmaker)
# BiocManager::install("DESeq2")
library(DESeq2)


# Path defintiions
# gtf <- "/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_10000_down_10000_longest_ends_05_09_2026/HUVEC/HUVEC_LR_SR_support_filtered.gtf"
# salmon_dir <- "/projects/splitorfs/work/short_RNA_seq_analysis/short_RNA_April_2025/HUVEC_TAMA_08_10_26/Salmon"
# gtf_name <- tools::file_path_sans_ext(basename(gtf))
# conditions <- "50NMD,5NMD,DHYPO,DMSO,DNOR"
# batch <- "2,3,4"

args <- commandArgs(trailingOnly = TRUE)
gtf <- args[1]
salmon_dir <- args[2]
conditions <- args[3]
batch <- args[4]
treat <- args[5]
control <- args[6]


gtf_name <- tools::file_path_sans_ext(basename(gtf))
deseq_outdir <- file.path(salmon_dir, "DESeq2")
dir.create(deseq_outdir, recursive = TRUE, showWarnings = FALSE)


get_name <- function(filename){
  new_name <- paste("Gibbs", filename, sep = "/")
  return(new_name)
}

create_filename <- function(filename){
 new_name <- paste(str_split_1(filename, "_")[1], str_split_1(filename, "_")[2], "merged", "fastp", "Gibbs", sep = "_")
  return(new_name)
}


setwd(salmon_dir)
sample_id <- lapply(list.files("Gibbs"), create_filename)
sample_id <- unlist(sample_id)
sample_id


conditions <- str_split(conditions, ",")
batch <- str_split(batch, ",")
batch <- rep(batch, length(conditions[[1]]))
batch <- unlist(batch)

files <- lapply(list.files(file.path(salmon_dir, "Gibbs")), get_name)
files <- unlist(files)
names(files) <- sample_id
conditions <- lapply(conditions, function(x) rep(x, each = 3))
conditions <- unlist(conditions)
conditions



txdb.filename <- paste(salmon_dir, paste0(gtf_name, ".sqlite"), sep = "/")
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)

txdb <- loadDb(txdb.filename)
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

colData <- DataFrame(
  conditions = conditions,
  batch     = batch
)
rownames(colData) <- sample_id



files_txi <- paste(files, "quant.sf", sep = "/")
txi <- tximport(files_txi, type="salmon", txOut=TRUE)


tx2gene <- data.frame(transcript_id = txdf[2], gene_id = txdf[1])
txi_g <- summarizeToGene(txi, tx2gene)

dds_gene <- DESeqDataSetFromTximport(txi_g, colData = colData, design = ~conditions+batch)



rld <- rlogTransformation(dds_gene)
pca <- plotPCA(rld, intgroup= c("conditions"))
svg(file.path(deseq_outdir, "pca_gene_plot.svg"), width = 7, height = 7)
plot(pca) 
dev.off()




# select conditions of interest
dds_gene_hypo <- dds_gene[, dds_gene$conditions %in% c(treat, control)]
# filter for expression in 4 samples and at least ten counts
dds_gene_hypo <- dds_gene_hypo[rowSums(counts(dds_gene_hypo)) >= 10, ]
keep <- rowSums(counts(dds_gene_hypo) >= 3) >= 3
dds_gene_hypo <- dds_gene_hypo[keep, ]
# reassign the reference level to mock
dds_gene_hypo$conditions <- droplevels(factor(dds_gene_hypo$conditions))
dds_gene_hypo$conditions <- relevel(dds_gene_hypo$conditions,
    ref = control
)
# rerun DeSEQ
dds_gene_hypo <- DESeq(dds_gene_hypo)

# perform differential testing
res_treat <- results(dds_gene_hypo,
    contrast = c("conditions", treat, control)
)
print(summary(res_treat))
# kept the independent filter, filter out NA values
res_treat <- res_treat[res_treat$baseMean > 0, ]
res_treat <- res_treat[complete.cases(res_treat), ]

# write output
write.csv(res_treat,
file=file.path(deseq_outdir, paste(treat, "vs", control, "DESeq2_results.csv", sep = "_")),
row.names=TRUE, quote=FALSE)

res_treat_padj <- res_treat[res_treat$padj < 0.05, ]
res_treat_padj <- na.omit(res_treat_padj)

write.csv(res_treat_padj,
file=file.path(deseq_outdir, paste(treat, "vs", control, "DESeq2_all_results_padj_smaller_0_05.csv", sep = "_")),
row.names=TRUE, quote=FALSE)

# filter by LFC as well:
res_treat_padj_lfc <- res_treat_padj[abs(res_treat_padj$log2FoldChange) > 0.58, ]

write.csv(res_treat_padj_lfc,
file=file.path(deseq_outdir, paste(treat, "vs", control, "DESeq2_all_results_padj_smaller_0_05_abs_LFC_greater_0_58.csv", sep = "_")),
row.names=TRUE, quote=FALSE)

# filter by LFC as well up and downregulated genes:
res_treat_padj_lfc_up <- res_treat_padj[res_treat_padj$log2FoldChange > 0.58, ]

write.csv(res_treat_padj_lfc_up,
file=file.path(deseq_outdir, paste(treat, "vs", control, "DESeq2_up_results_padj_smaller_0_05_abs_LFC_greater_0_58.csv", sep = "_")),
row.names=TRUE, quote=FALSE)


# filter by LFC as well up and downregulated genes:
res_treat_padj_lfc_down <- res_treat_padj[res_treat_padj$log2FoldChange < -0.58, ]

write.csv(res_treat_padj_lfc_down,
file=file.path(deseq_outdir, paste(treat, "vs", control, "DESeq2_down_results_padj_smaller_0_05_abs_LFC_greater_0_58.csv", sep = "_")),
row.names=TRUE, quote=FALSE)

utils::sessionInfo()