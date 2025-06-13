library(DESeq2)
library(stringr)
library(tximportData)
library(GenomicFeatures)

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]


################################################################################
# READ IN QUANT.sf files                                                       #
################################################################################

gtf <- "/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.gtf"
txdb.filename <- "/projects/splitorfs/work/reference_files/clean_Ensembl_ref/Ensembl_equality_and_TSL_filtered.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)

txdb <- loadDb(txdb.filename)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")



get_name <- function(filename){
  new_name <- file.path(in_dir, filename, "quant.sf")
  return(new_name)
}

files <- lapply(list.files(in_dir), get_name)
files <- unlist(files)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

idx_files <- c()
if (!(identical(
    list.files(in_dir,
        pattern = paste0("*", "idxstats.out")
    ),
    character(0)
))) {
    idx_files <- c(
        idx_files,
        paste0(
            in_dir,
            "/",
            list.files(in_dir, pattern = paste0("*", "idxstats.out"))
        )
    )
}

merged_counts_df <- NULL
for (file in idx_files) {
    sample <- basename(file)
    sample <- substr(sample, 20, 43)
    sample <- str_split(sample, "\\.")[[1]][1]

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


################################################################################
# FILTER FOR mRNA mapping reads                                                #
################################################################################
MANE_mRNAs <- read.table(
    "/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_rna.fna.fai",
    header = FALSE,
    sep = "\t",
    col.names = c("tID", "length", "some_nr", "some_other_nr", "yet_another_nr")
)

merged_counts_df_mRNA <- merged_counts_df[
    rownames(merged_counts_df) %in% MANE_mRNAs$tID,
]

merged_counts_matrix_mRNA <- data.matrix(merged_counts_df_mRNA)
merged_counts_matrix_mRNA <- merged_counts_matrix_mRNA[
    rowSums(merged_counts_matrix_mRNA) > 10,
]

col_data <- data.frame(
    row.names = colnames(merged_counts_df_mRNA),
    condition = factor(c(
        "SeRP_CHX_A2",
        "SeRP_CHX_A2",
        "SeRP_CHX_A2",
        "SeRP_CHX_A2",
        "SeRP_Puro_A2",
        "SeRP_Puro_A2",
        "SeRP_Puro_A2",
        "SeRP_CHX_B1",
        "SeRP_CHX_B1",
        "SeRP_CHX_B1",
        "SeRP_CHX_B1",
        "SeRP_Puro_B1",
        "SeRP_Puro_B1",
        "SeRP_Puro_B1",
        "In_CHX",
        "In_CHX",
        "In_CHX",
        "In_CHX",
        "In_Puro",
        "In_Puro",
        "In_Puro",
        "SeRP_CHX_M",
        "SeRP_CHX_M",
        "SeRP_CHX_M",
        "SeRP_CHX_M"
    )),
    batch = factor(c(
    "CHX_1", "CHX_2", "CHX_3", "CHX_4",
    "Puro_1", "Puro_2", "Pruo_3",
    "CHX_1", "CHX_2", "CHX_3", "CHX_4",
    "Puro_1", "Puro_2", "Pruo_3",
    "CHX_1", "CHX_2", "CHX_3", "CHX_4",
    "Puro_1", "Puro_2", "Pruo_3",
    "CHX_1", "CHX_2", "CHX_3", "CHX_4"
    ))
)



################################################################################
# Variance transformation and plotting                                         #
################################################################################

dds_mRNA <- DESeqDataSetFromMatrix(merged_counts_df_mRNA,
    colData = col_data,
    design = ~condition
)
vid_mRNA <- vst(dds_mRNA, blind = FALSE)

png(
    file.path(
        in_dir,
        "PCA_transcriptome_counts_greater_10_vid_mRNA.png"
    ),
    width = 800,
    height = 800
)
print(plotPCA(vid_mRNA, intgroup = c("condition")))
dev.off()


# RLOG TRANSFORMATION
rid_mRNA <- rlog(dds_mRNA, blind = FALSE)
print("correlation rid")
cor(assay(rid_mRNA))

png(
    file.path(
        in_dir,
        "PCA_transcriptome_counts_greater_10_rlog_mRNA.png"
    ),
    width = 800,
    height = 800
)
print(plotPCA(rid_mRNA, intgroup = c("condition")))
dev.off()



# MA plot
dds_mRNA <- DESeq(dds_mRNA)
res_mRNA <- results(dds_mRNA)

png(file.path(in_dir, "MA_plot_transcriptome_mRNA.png"),
    width = 800, height = 800
)
print(plotMA(res_mRNA), ylim = c(-2, 2))
dev.off()


resLFC_mRNA <- lfcShrink(dds_mRNA, coef = 2, type = "apeglm")
png(file.path(in_dir, "MA_res_LFC_plot_transcriptome_mRNA.png"),
    width = 800, height = 800
)
print(plotMA(resLFC_mRNA, ylim = c(-2, 2)))
dev.off()


################################################################################
# OBTAINING DEGs (transcripts) between SeRP and Input ##########################
################################################################################
perform_diff_analysis <- function(dds_mRNA, treatment, control, outname, in_dir)
 {
    # select conditions of interest
    dds_mRNA_treat <- dds_mRNA[, dds_mRNA$condition %in% c(treatment, control)]
    # filter for expression in 4 samples and at least ten counts
    dds_mRNA_treat <- dds_mRNA_treat[rowSums(counts(dds_mRNA_treat)) >= 10, ]
    keep <- rowSums(counts(dds_mRNA_treat) >= 1) >= 4
    dds_mRNA_treat <- dds_mRNA_treat[keep, ]
    # reassign the reference level to mock
    dds_mRNA_treat$condition <- droplevels(dds_mRNA_treat$condition)
    dds_mRNA_treat$condition <- relevel(dds_mRNA_treat$condition,
        ref = control
    )
    # rerun DeSEQ
    dds_mRNA_treat <- DESeq(dds_mRNA_treat)

    # perform differential testing
    res_treat <- results(dds_mRNA_treat,
        contrast = c("condition", treatment, control)
    )
    # kept the independent filter, filter out NA values
    res_treat <- res_treat[res_treat$baseMean > 0, ]
    res_treat_pval <- res_treat[res_treat$pvalue < 0.05, ]
    res_treat_pval <- na.omit(res_treat_pval)
    # write output
    writeLines(
        rownames(res_treat_pval[
            (res_treat_pval$padj < 0.05) & (res_treat_pval$log2FoldChange > 0.5),
        ]),
        file.path(in_dir, "DEGs", outname)
    )

    writeLines(
        rownames(res_treat_pval[
            (res_treat_pval$padj < 0.05) & (res_treat_pval$log2FoldChange < -0.5),
        ]),
        file.path(in_dir, "DEGs", paste(str_split(outname, "\\.")[[1]][1], "downreg.txt", sep = "_"))
    )

    svg(file.path(in_dir, "DEGs", paste(treatment, control, "MA_plot.svg", sep ="_")))
    print(plotMA(res_treat), ylim = c(-2, 2))
    dev.off()

    png(file.path(in_dir, "DEGs", paste(treatment, control, "MA_plot.png", sep ="_")))
    print(plotMA(res_treat), ylim = c(-2, 2))
    dev.off()


    return(res_treat_pval)
}

perform_diff_analysis(dds_mRNA, "SeRP_CHX_A2", "SeRP_CHX_M", "DEGs_CHX_A2_mock_mRNA.txt", in_dir)
perform_diff_analysis(dds_mRNA, "SeRP_CHX_B1", "SeRP_CHX_M", "DEGs_CHX_B1_mock_mRNA.txt", in_dir)
perform_diff_analysis(dds_mRNA, "SeRP_Puro_A2", "In_Puro", "DEGs_Puro_A2_mRNA.txt", in_dir)
perform_diff_analysis(dds_mRNA, "SeRP_Puro_B1", "In_Puro", "DEGs_Puro_B1_mRNA.txt", in_dir)
perform_diff_analysis(dds_mRNA, "SeRP_CHX_M", "In_CHX", "DEGs_Mock_In_CHX_mRNA.txt", in_dir)
perform_diff_analysis(dds_mRNA, "SeRP_CHX_A2", "In_CHX", "DEGs_A2_In_CHX_mRNA.txt", in_dir)
perform_diff_analysis(dds_mRNA, "SeRP_CHX_B1", "In_CHX", "DEGs_B1_In_CHX_mRNA.txt", in_dir)


