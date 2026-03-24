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


# args <- commandArgs(trailingOnly = TRUE)
nmd_path <- "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD"
ri_path <- "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI"
coarse <- TRUE



cell_type_to_sample_dict <- list(
    "CM_IPSC" = c("ERR3367797", "ERR3367798"),
    "HUVEC_untreated" = c("OHMX20220060_001", "OHMX20220060_002", "OHMX20220060_003"),
    "HUVEC_treated" = c("OHMX20220060_004", "OHMX20220060_005", "OHMX20220060_006"),
    "HEK control" = c("HCT_N_AID_UPF1_0h_IAA_1Ens_110",
                        "HCT_N_AID_UPF1_0h_IAA_2Ens_110", 
                        "HCT_N_AID_UPF1_0h_IAA_3Ens_110"),
    "healthy breast hmec" =  c("SRR8590754", 
                                "SRR8590755", 
                                "SRR8590756",
                                "SRR8590757",
                                "SRR8590758", 
                                "SRR8590759"),
    "MCF 10A fibrocystic" = c("SRR8590766",
                                "SRR8590767", 
                                "SRR8590768",
                                "SRR8590769"),
    "T47D metastatic" = c("SRR8590778", 
                            "SRR8590779", 
                            "SRR8590780", 
                            "SRR8590781", 
                            "SRR8590782", 
                            "SRR8590783"),
    "ZR75-1 metastatic" = c("SRR8590790"),
    "glioblastoma U251" = c("SRR10533453", 
                            "SRR10533454", 
                            "SRR10533455", 
                            "SRR10533456",
                            "SRR10533457", 
                            "SRR10533458", 
                            "SRR10533459", 
                            "SRR105334560"),
     "glioblastoma U343" = c("SRR10533461",
                            "SRR10533462", 
                            "SRR10533463", 
                            "SRR10533464", 
                            "SRR10533465", 
                            "SRR10533466",
                            "SRR10533467", 
                            "SRR10533468", 
                            "SRR10533469"),
    "HUVEC_normoxia" = c("uf_muellermcnicoll_2025_04_01_huvec_dnor_2", "uf_muellermcnicoll_2025_04_02_huvec_dnor_3", "uf_muellermcnicoll_2025_04_03_huvec_dnor_4"),
    "HUVEC_hypoxia" = c("uf_muellermcnicoll_2025_04_04_huvec_dhypo_2", "uf_muellermcnicoll_2025_04_05_huvec_dhypo_3","uf_muellermcnicoll_2025_04_06_huvec_dhypo_4")
)



cell_type_to_sample_dict_coarse <- list(
    "CM_IPSC" = c("ERR3367797", "ERR3367798"),
    "HUVEC" = c("OHMX20220060_001", "OHMX20220060_002", "OHMX20220060_003", 
                "OHMX20220060_004", "OHMX20220060_005", "OHMX20220060_006",
                "uf_muellermcnicoll_2025_04_01_huvec_dnor_2",
                "uf_muellermcnicoll_2025_04_02_huvec_dnor_3",
                "uf_muellermcnicoll_2025_04_03_huvec_dnor_4", 
                "uf_muellermcnicoll_2025_04_04_huvec_dhypo_2", 
                "uf_muellermcnicoll_2025_04_05_huvec_dhypo_3",
                "uf_muellermcnicoll_2025_04_06_huvec_dhypo_4"),
    "HEK control" = c("HCT_N_AID_UPF1_0h_IAA_1Ens_110", "HCT_N_AID_UPF1_0h_IAA_2Ens_110", "HCT_N_AID_UPF1_0h_IAA_3Ens_110"),
    "breast cancer" = c("SRR8590754", "SRR8590755", "SRR8590756", "SRR8590757",
                "SRR8590758", "SRR8590759", "SRR8590766", "SRR8590767", "SRR8590768", "SRR8590769",
                "SRR8590778", "SRR8590779", "SRR8590780", "SRR8590781", "SRR8590782", "SRR8590783",
        "SRR8590790"),
    "glioblastoma" = c("SRR10533453", "SRR10533454", "SRR10533455", "SRR10533456",
                "SRR10533457", "SRR10533458", "SRR10533459", "SRR105334560", "SRR105334561",
                "SRR105334562", "SRR105334563", "SRR105334564", "SRR105334565", "SRR105334566",
                "SRR105334567", "SRR105334568", "SRR105334569")
)


#################################################################################
# ------------------ HELPER FUNCTION DEFINITIONS             ------------------ #
#################################################################################
#################################################################################
# ------------------ HELPER FUNCTION DEFINITIONS             ------------------ #
#################################################################################
upset_list_gen_regions_per_cell_type <- function(cell_type_to_sample_dict, validated_so_df){
    cell_type_concat_dfs <- list()
    upsetlist_genomic <- list()
    for (cell_type in names(cell_type_to_sample_dict)) {
        samples_of_cell_type <- cell_type_to_sample_dict[cell_type]
        samples_of_cell_type <- unlist(samples_of_cell_type)[unlist(samples_of_cell_type) %in% colnames(validated_so_df)]
        keep_cols <- c(samples_of_cell_type, c("ID"))
        df_sample <- validated_so_df[,keep_cols]

        if (is.data.frame(df_sample) && nrow(df_sample) > 0){
            df_sample_validated <- df_sample[rowSums(df_sample == "True", na.rm = TRUE) > 0, "ID"]
            if (length(df_sample_validated > 0)) {
            upsetlist_genomic[[cell_type]] <- unique(df_sample_validated) 
            }
        }
    }
    return(upsetlist_genomic)
}


get_upset_list_from_csv <- function(path, coarse) {
    file_path <-  file.path(path, "val_dna_overlapping_ur_df.csv")
    validated_so_df <- read.csv(file_path, sep = ",")

    if (coarse){
        upsetlist_genomic <- upset_list_gen_regions_per_cell_type(cell_type_to_sample_dict_coarse, validated_so_df)
    } else {
        upsetlist_genomic <- upset_list_gen_regions_per_cell_type(cell_type_to_sample_dict, validated_so_df)
    }
    return(upsetlist_genomic)
}

upsetlist_genomic_nmd <- get_upset_list_from_csv(nmd_path, coarse)
upsetlist_genomic_ri <- get_upset_list_from_csv(ri_path, coarse)

sorted_samples_nmd <- sort(names(upsetlist_genomic_nmd))
sorted_samples_ri <- sort(names(upsetlist_genomic_ri))


if (coarse) {
    nsets <- 4
    set_size.scale_max_nmd <- 150
    set_size.scale_max_ri <- 300
    nmd_out <- 'upset_NMD_noNMD_24_03_26_coarse.svg'
    ri_out <-'upset_RI_noNMD_24_03_26_coarse.svg'
} else {
    nsets <- 7
    set_size.scale_max_nmd <- 150
    set_size.scale_max_ri <- 300
    nmd_out <- 'upset_NMD_noNMD_24_03_26.svg'
    ri_out <-'upset_RI_noNMD_24_03_26.svg'
}

svg(file.path(nmd_path, nmd_out), width = 14, height = 6)
    upset(fromList(upsetlist_genomic_nmd), order.by = "freq", nsets = nsets, point.size = 2.25, line.size = 1.5, set_size.show = TRUE,
        set_size.scale_max = set_size.scale_max_nmd, nintersects=30, 
        mainbar.y.label = "NMD Unique Region Intersections", sets.x.label = "Matching NMD unique regions per cell type",
        text.scale = 1.7)
dev.off()
svg(file.path(ri_path, ri_out), width = 14, height = 6)
    upset(fromList(upsetlist_genomic_ri), order.by = "freq", nsets = nsets, point.size = 2.25, line.size = 1.5, set_size.show = TRUE,
        set_size.scale_max = set_size.scale_max_ri, sets = sorted_samples_ri, nintersects = 30,
        mainbar.y.label = "RI Unique Region Intersections", sets.x.label = "Matching RI unique regions per cell type",
        text.scale = 1.7)
dev.off()

utils::sessionInfo()
