if (!requireNamespace("types", quietly = TRUE)) {
    install.packages("types", repos = "https://cloud.r-project.org/")
  }

if (!requireNamespace("UpSetR", quietly = TRUE)) {
    install.packages("UpSetR", repos = "https://cloud.r-project.org/")
  }

if (!requireNamespace("lemon", quietly = TRUE)) {
    install.packages("lemon", repos = "https://cloud.r-project.org/")
  }

if (!requireNamespace("cmapR", quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")

    BiocManager::install("cmapR")
  }
if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}
if (!requireNamespace("readr", quietly = TRUE)) {
    install.packages("readr")
}


library(glue)
library(types)
library(stringr)
library(dplyr)
library(UpSetR)
library(knitr)
library(lemon)
library(cmapR)
library(data.table)


library(readr)
path_NMD <- '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample/NMD_genome'

samples_to_display <- c('OHMX20220060_001_genomicunique_regions.csv',
 'OHMX20220060_003_genomicunique_regions.csv',
 'OHMX20220060_004_genomicunique_regions.csv',
 'OHMX20220060_006_genomicunique_regions.csv',
 'ERR3367797_genomicunique_regions.csv',
 'ERR3367798_genomicunique_regions.csv')

upsetlist_genomic <- list()
for (unique_regions in samples_to_display) {
    data <- read_csv(file.path(path_NMD, unique_regions), col_names = TRUE) %>% select(-1)
    data <- data[data$significant == 1,]
    upsetlist_genomic <- c(upsetlist_genomic, list(data$ID))
    }

names(upsetlist_genomic) <- c('HUVEC_untreated_1', 'HUVEC_untreated_2', 'HUVEC_treated_1', 'HUVEC_treated_2', 'CM_IPSC_1', 'CM_IPSC_2')

svg(file.path(path_NMD, 'upset_NMD_poster_31_03_25.svg'), width = 12, height = 6)
upsetlist_genomic <- lapply(upsetlist_genomic, function(x) x)
upset(fromList(upsetlist_genomic), order.by = "freq", nsets = 12, point.size = 2.25, line.size = 1.5, set_size.show = TRUE,
      set_size.scale_max = 350, 
      mainbar.y.label = "Unique Region Intersections", sets.x.label = "Matching unique regions per dataset",
      text.scale = c(1, 1, 0.9, 1, 1, 1))
dev.off()
