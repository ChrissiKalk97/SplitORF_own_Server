library(DESeq2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]


################################################################################
# READ IN IDX STATS                                                            #
################################################################################
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
    "/projects/serp/work/references/Chaetomium_thermophilum_protein_coding.fasta.fai",
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
                            "In_WT",
                            "In_WT",
                            "In_WT",
                            "In_WT",

                            "In_S",
                            "In_S",
                            "In_S",
                            "In_S",

                            "E_WT",
                            "E_WT",
                            "E_WT",
                            "E_WT",

                            "E_S",
                            "E_S",
                            "E_S",
                            "E_S",


                            "2025_05_17_MD_In_W1_ND"
    )),
    batch = factor(c(
    "WT_1", "WT_2", "WT_3", "WT_4",
    "S_1", "S_2", "S_3", "S_4",
    "WT_1", "WT_2", "WT_3", "WT_4",
    "S_1", "S_2", "S_3", "S_4",
    "WT_1"
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
perform_diff_analysis <- function(dds_mRNA, treatment, control, outname, in_dir, lfc = 0.5, p_adj = 0.05)
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
    res_treat <- res_treat[complete.cases(res_treat), ]
    res_treat_pval <- res_treat[res_treat$pvalue < p_adj, ]
    res_treat_pval <- na.omit(res_treat_pval)
    # write output
    writeLines(
        rownames(res_treat_pval[
            (res_treat_pval$padj < p_adj) & (res_treat_pval$log2FoldChange > lfc),
        ]),
        file.path(in_dir, "DEGs", paste(str_split(outname, "\\.")[[1]][1], paste(p_adj, "txt", sep = "."), sep = "_")
        )
    )


    write.csv(res_treat_pval,
    file=file.path(in_dir, "DEGs", paste(str_split(outname, "\\.")[[1]][1], "all_results_p_smaller", paste(p_adj, "csv", sep = "."), sep = "_")),
    row.names=TRUE, quote=FALSE)

    writeLines(
        rownames(res_treat_pval[
            (res_treat_pval$padj < p_adj) & (res_treat_pval$log2FoldChange < -lfc),
        ]),
        file.path(in_dir, "DEGs", paste(str_split(outname, "\\.")[[1]][1], p_adj, "downreg.txt", sep = "_"))
    )

    svg(file.path(in_dir, "DEGs", paste(treatment, control, p_adj, "MA_plot.svg", sep ="_")))
    print(plotMA(res_treat), ylim = c(-2, 2))
    dev.off()

    png(file.path(in_dir, "DEGs",  paste(treatment, control, p_adj, "MA_plot.png", sep ="_")))
    print(plotMA(res_treat), ylim = c(-2, 2))
    dev.off()


    return(res_treat_pval)
}


perform_diff_analysis(dds_mRNA, "E_WT", "In_WT", paste0("E_over_In_WT", "_0_5",".txt"), in_dir, 0.5)
perform_diff_analysis(dds_mRNA, "E_S", "In_S", paste0("E_over_In_S", "_0_5",".txt"), in_dir, 0.5)
perform_diff_analysis(dds_mRNA, "E_S", "E_WT", paste0("E_S_over_E_WT", "_0_5",".txt"), in_dir, 0.5)


perform_diff_analysis(dds_mRNA, "E_WT", "In_WT", paste0("E_over_In_WT", "_1_0",".txt"), in_dir, 1.0, 0.01)
perform_diff_analysis(dds_mRNA, "E_S", "In_S", paste0("E_over_In_S", "_1_0",".txt"), in_dir, 1.0, 0.01)
perform_diff_analysis(dds_mRNA, "E_S", "E_WT", paste0("E_S_over_E_WT", "_1_0",".txt"), in_dir, 1.0, 0.01)