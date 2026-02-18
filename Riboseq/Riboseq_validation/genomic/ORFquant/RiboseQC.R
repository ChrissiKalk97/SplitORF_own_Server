library(RiboseQC)


args <- commandArgs(trailingOnly = TRUE)
twobit <- args[1]
gtf <- args[2]
outdir <- args[3]
genome_fasta <- args[4]
bam_path <- args[5]


# file.exists(twobit)
# prepare_annotation_files(annotation_directory = outdir,
#                          twobit_file = twobit,
#                          gtf_file = gtf,
#                          scientific_name = "Homo.sapiens",
#                          annotation_name = "Ensembl110",
#                          export_bed_tables_TxDb = F,
#                          forge_BSgenome = F,
#                          create_TxDb = T)

load_annotation( paste(outdir, "Homo_sapiens.GRCh38.110.chr.no.comment.gtf_Rannot", sep = "/"))

Riboseq_bams <- list.files(path = bam_path, pattern = "\\.bam$", full.names = TRUE)

print(Riboseq_bams)

sample_names <- character(length(Riboseq_bams))  # preallocate for speed

for (i in seq_along(Riboseq_bams)) {
    bam <- Riboseq_bams[i]
    sample <- basename(bam)
    if (grepl("NMD_sorted", sample)) {
        sample <- sub("\\_NMD_sorted\\.bam$", "", sample)
    } else if (grepl("_dedup_filtered_sorted", sample)) {
        sample <- sub("\\_dedup\\_filtered\\_sorted.bam", "", sample)
    }
    sample_names[i] <- file.path(outdir, "RiboseQC", "prepared_bam_files", sample)
}

print(length(sample_names))
print(sample_names)

print(length(Riboseq_bams))

i = 1
for (bam_file in Riboseq_bams) {
    sample <- basename(bam_file)
    sample <- sub("\\_NMD_sorted\\.bam$", "", sample)
    RiboseQC_analysis(annotation_file=paste(outdir, "Homo_sapiens.GRCh38.110.chr.no.comment.gtf_Rannot", sep = "/"),
    bam_files = bam_file,
    report_file = file.path(outdir, "RiboseQC", paste(sample, "04_08_25.html", sep = "_")),
    write_tmp_files = F,
    # sample_names = sample_names[i],
    dest_names = sample_names[i]
    )
    i = i + 1
}


