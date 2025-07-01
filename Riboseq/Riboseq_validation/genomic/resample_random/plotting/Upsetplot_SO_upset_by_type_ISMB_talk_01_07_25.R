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


args <- commandArgs(trailingOnly = TRUE)
nmd_path <- args[1]
ri_path <- args[2]


sample_to_cell_type_dict <- c(
    "ERR3367797" = "CM_IPSC",
    "ERR3367798" = "CM_IPSC",
    "OHMX20220060_001" = "HUVEC_untreated",
    "OHMX20220060_002" = "HUVEC_untreated",
    "OHMX20220060_003" = "HUVEC_untreated",
    "OHMX20220060_004" = "HUVEC_treated",
    "OHMX20220060_005" = "HUVEC_treated",
    "OHMX20220060_006" = "HUVEC_treated",
    "SRR11294608" = "Leukemia_untreated",
    "SRR11294609" = "Leukemia_untreated",
    "SRR11294610"  = "Leukemia_treated",
    "SRR11294611"  = "Leukemia_treated",
    "uf_muellermcnicoll_2025_04_01_huvec_dnor_2" = "HUVEC_normoxia",
    "uf_muellermcnicoll_2025_04_02_huvec_dnor_3" = "HUVEC_normoxia",
    "uf_muellermcnicoll_2025_04_03_huvec_dnor_4" = "HUVEC_normoxia",
    "uf_muellermcnicoll_2025_04_04_huvec_dhypo_2" = "HUVEC_hypoxia",
    "uf_muellermcnicoll_2025_04_05_huvec_dhypo_3" = "HUVEC_hypoxia",
    "uf_muellermcnicoll_2025_04_06_huvec_dhypo_4" = "HUVEC_hypoxia"
)

cell_type_to_sample_dict <- list(
    "CM_IPSC" = c("ERR3367797", "ERR3367798"),
    "HUVEC_untreated" = c("OHMX20220060_001", "OHMX20220060_002", "OHMX20220060_003"),
    "HUVEC_treated" = c("OHMX20220060_004", "OHMX20220060_005", "OHMX20220060_006"),
    "Leukemia_untreated" = c("SRR11294608", "SRR11294609"),
    "Leukemia_treated" = c("SRR11294610", "SRR11294611"),
    "HUVEC_normoxia" = c("uf_muellermcnicoll_2025_04_01_huvec_dnor_2", "uf_muellermcnicoll_2025_04_02_huvec_dnor_3", "uf_muellermcnicoll_2025_04_03_huvec_dnor_4"),
    "HUVEC_hypoxia" = c("uf_muellermcnicoll_2025_04_04_huvec_dhypo_2", "uf_muellermcnicoll_2025_04_05_huvec_dhypo_3","uf_muellermcnicoll_2025_04_06_huvec_dhypo_4")
)


#################################################################################
# ------------------ HELPER FUNCTION DEFINITIONS             ------------------ #
#################################################################################

get_validated_urs_df <- function(validated_regions_df){
    # remove the index column from the dataframes
    rownames(validated_regions_df) <- validated_regions_df[[1]]
    validated_regions_df[[1]] <- NULL

    # subset for significant regions only
    validated_regions_df <- validated_regions_df[validated_regions_df$significant == 1]
    validated_regions_df <- validated_regions_df %>% mutate(genomic_region = paste(chr_unique, start, stop, sep = "_"))

    return(validated_regions_df)
}


upset_list_gen_regions_per_cell_type <- function(cell_type_to_sample_dict, validated_dfs){
    cell_type_concat_dfs <- list()
    upsetlist_genomic <- list()
    for (cell_type in names(cell_type_to_sample_dict)) {
        samples_of_cell_type <- cell_type_to_sample_dict[cell_type]
        dfs_to_combine <- validated_dfs[unlist(samples_of_cell_type)]
        combined_df <- bind_rows(dfs_to_combine, .id = " ID")
        print(length(rownames(combined_df)))
        cell_type_concat_dfs[[cell_type]] <- combined_df
        upsetlist_genomic[[cell_type]] <- unique(combined_df$genomic_region)
        # print(unique(combined_df$genomic_region))
        print(length(unique(combined_df$genomic_region)))
    }
    return(upsetlist_genomic)
}


get_upset_list_from_csvs <- function(path) {
    files <- list.files(path, pattern = 
                "*_genomicunique_regions.csv"
            )
    file_paths <-  file.path(path, files)
    so_with_riboseq_dfs <- lapply(file_paths,
            fread,
            header = TRUE,
            sep = ","
        )

    validated_dfs <- lapply(so_with_riboseq_dfs, get_validated_urs_df)

    sample_names <- sub("_genomicunique_regions\\.csv$", "", files)

    names(validated_dfs) <- sample_names

    upsetlist_genomic <- upset_list_gen_regions_per_cell_type(cell_type_to_sample_dict, validated_dfs)
    return(upsetlist_genomic)

}

upsetlist_genomic_nmd <- get_upset_list_from_csvs(nmd_path)
upsetlist_genomic_ri <- get_upset_list_from_csvs(ri_path)

sorted_samples_nmd <- sort(names(upsetlist_genomic_nmd))
sorted_samples_ri <- sort(names(upsetlist_genomic_ri))


svg(file.path(nmd_path, 'upset_NMD_cell_types_01_07_25.svg'), width = 12, height = 6)
upset(fromList(upsetlist_genomic_nmd), order.by = "freq", nsets = 7, point.size = 2.25, line.size = 1.5, set_size.show = TRUE,
      set_size.scale_max = 500,
      mainbar.y.label = "NMD Unique Region Intersections", sets.x.label = "Matching NMD unique regions per cell type",
      text.scale = c(1, 1, 0.9, 1, 1, 1))
dev.off()


svg(file.path(ri_path, 'upset_RI_cell_types_01_07_25.svg'), width = 12, height = 6)
upset(fromList(upsetlist_genomic_ri), order.by = "freq", nsets = 7, point.size = 2.25, line.size = 1.5, set_size.show = TRUE,
      set_size.scale_max = 300, sets = sorted_samples_ri,
      mainbar.y.label = "RI Unique Region Intersections", sets.x.label = "Matching RI unique regions per cell type",
      text.scale = c(1, 1, 0.9, 1, 1, 1))
dev.off()
