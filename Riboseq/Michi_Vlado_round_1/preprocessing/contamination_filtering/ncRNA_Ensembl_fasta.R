BiocManager::install("biomaRt")
library(biomaRt)
listEnsemblArchives()

ensembl_110 <- useEnsembl(biomart = 'genes', 
                       dataset = 'hsapiens_gene_ensembl',
                       version = 110,
                       mirror = "www")

# mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

biotypes <- c("misc_RNA", "rRNA", "ribozym", "sRNA", "scRNA", "scaRNA", "snRNA", "snoRNA", "vaultRNA")

data <- getBM(attributes=c("ensembl_transcript_id", "transcript_biotype"),
              filters = "transcript_biotype",
              values = biotypes,
              mart = ensembl_110)

nc_sequences <- getSequence(type = "ensembl_transcript_id",
              seqType = "cdna",
              id = data[[1]],
              mart = ensembl_110)



exportFASTA(nc_sequences, "/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/ncRNA/redownload/Ensembl_all_nc_biotypes.fasta")