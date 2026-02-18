library(ORFquant)
library(BSgenome)

args <- commandArgs(trailingOnly = TRUE)
twobit <- args[1]
gtf <- args[2]
outdir <- args[3]
genome_fasta <- args[4]
for_orfquant_path <- args[5]

print(typeof(twobit))
print(length(twobit))
print(gtf)
print(outdir)

# trace(".copyTwobitFile", where = asNamespace("BSgenome"), tracer = quote(print(list(twobit, file.path(outdir, "BSgenome.Homo.sapiens.Ensembl110", "inst", "extdata"), outdir))))

# file.exists(twobit)
# prepare_annotation_files(annotation_directory = outdir,
#                          twobit_file = twobit,
#                          gtf_file = gtf,
#                          scientific_name = "Homo.sapiens",
#                          annotation_name = "Ensembl110",
#                          genome_seq = genome_fasta,
#                          export_bed_tables_TxDb = F,
#                          forge_BSgenome = F,
#                          create_TxDb = T)


annotation_Rannot <- file.path(outdir, "Homo_sapiens.GRCh38.110.chr.no.comment.gtf_Rannot")
load_annotation(annotation_Rannot )


Riboseq_orfquant_files <- list.files(path = for_orfquant_path, pattern = "_for_ORFquant$", full.names = TRUE)
print(Riboseq_orfquant_files)
for (file in Riboseq_orfquant_files)
    {
        print(file)
        sample <- basename(file)
        sample <- sub("_for_ORFquant$", "", sample)
        result_dir <- file.path(outdir, "ORFquant_output", sample)
        print(result_dir)
        run_ORFquant(for_ORFquant_file = file,
           annotation_file = annotation_Rannot ,
           n_cores = 16,
           prefix = result_dir)
    }
