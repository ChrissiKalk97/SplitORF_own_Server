library(ORFquant)
library(BSgenome)

args <- commandArgs(trailingOnly = TRUE)
twobit <- args[1]
gtf <- args[2]
outdir <- args[3]
genome_fasta <- args[4]
# bam_path <- args[4]

print(typeof(twobit))
print(length(twobit))
print(gtf)
print(outdir)

# trace(".copyTwobitFile", where = asNamespace("BSgenome"), tracer = quote(print(list(twobit, file.path(outdir, "BSgenome.Homo.sapiens.Ensembl110", "inst", "extdata"), outdir))))

file.exists(twobit)
prepare_annotation_files(annotation_directory = outdir,
                         twobit_file = twobit,
                         gtf_file = gtf,
                         scientific_name = "Homo.sapiens",
                         annotation_name = "Ensembl110",
                         genome_seq = genome_fasta,
                         export_bed_tables_TxDb = F,
                         forge_BSgenome = F,
                         create_TxDb = T)

