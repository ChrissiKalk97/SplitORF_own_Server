# install.packages("devtools")
# library(devtools)
# install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)

library(riboWaltz)
library(Biostrings) 
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]

plot_dir = file.path(in_dir, 'Ribowaltz')
gtf_file <- '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf'
annotation_df <- create_annotation(gtfpath = gtf_file, dataSource = "RefSeqMANE0.95", organism = "Homo sapiens")


setwd(plot_dir)

# Read in the FASTA file (DNA sequences example)
fasta_data <- readDNAStringSet("/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta")

# Retrieve the headers (sequence names)
headers <- names(fasta_data)


####################################################
#READ in BAM FILES of deduplicated reads
####################################################


reads_list <- bamtolist(bamfolder = "/projects/serp/work/Output/April_2025/importins/transcriptome_mapping/bowtie1/filtered/q10", annotation = annotation_df)

join_sample_name <- function(string_vec){
  string_vec[9] <- str_split(string_vec[9], '\\.')[[1]]
  return(paste(string_vec[3:9], collapse = '_'))
}

new_sample_names <- sapply(str_split(names(reads_list), '_'), join_sample_name)
names(reads_list) <- new_sample_names



psite_offset <- psite(reads_list, flanking = 6, extremity = "auto", plot = TRUE)

reads_psite_list <- psite_info(reads_list, psite_offset)


codon_coverage_example <- codon_coverage(reads_psite_list, annotation_df, psite = FALSE)


#######################################################################################################
# PLOTTING
#######################################################################################################

#### LENGTH DISTRIBUTION ##############################################################################
dir_path <- file.path(plot_dir, "read_lengths")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}


for (sample in names(reads_psite_list)) {
  sample_dedup_list <- list()
  sample_dedup_list[[sample]] <- c(sample)
  example_length_dist <- rlength_distr(reads_psite_list, 
                                       sample = sample_dedup_list,
                                       cl = 99, 
                                       colour = "#333f50")
  
  plot_name <- paste("plot", sample, sep = '_')
  png(file.path(plot_dir, "read_lengths", paste(sample, 'read_lengths.png', sep = '_')), width = 1000, height = 500)
  print(example_length_dist[[plot_name]])
  dev.off()
}


#define conditions and replicates
HEK_samples <- list('In_Puro' = c('2025_05_36_RR_In_Puro_1',
                                  '2025_05_37_RR_In_Puro_2',
                                  '2025_05_38_RR_In_Puro_3'),
                      'Mock_CHX' = c('2025_05_39_RR_M_CHX_E1',
                                    '2025_05_40_RR_M_CHX_E2',
                                    '2025_05_41_RR_M_CHX_E3',
                                    '2025_05_42_RR_M_CHX_E4'),
                      'In_CHX' = c('2025_05_32_RR_In_CHX_1',
                                        '2025_05_33_RR_In_CHX_2',
                                        '2025_05_34_RR_In_CHX_3',
                                        '2025_05_35_RR_In_CHX_4'),
                      'A2_CHX' = c('2025_05_18_RR_A2_CHX_E1',
                                    '2025_05_19_RR_A2_CHX_E2',
                                    '2025_05_20_RR_A2_CHX_E3',
                                    '2025_05_21_RR_A2_CHX_E4'),
                      'A2_Puro' = c('2025_05_22_RR_A2_Puro_E1',
                                    '2025_05_23_RR_A2_Puro_E2',
                                    '2025_05_24_RR_A2_Puro_E3'),
                      'B1_CHX' = c('2025_05_25_RR_B1_CHX_E1',
                                   '2025_05_26_RR_B1_CHX_E2',
                                   '2025_05_27_RR_B1_CHX_E3',
                                   '2025_05_28_RR_B1_CHX_E4'),
                      'B1_Puro' = c('2025_05_29_RR_B1_Puro_E1',
                                    '2025_05_30_RR_B1_Puro_E2',
                                    '2025_05_31_RR_B1_Puro_E3')
                      )


# read length distribution
example_length_dist <- rlength_distr(reads_psite_list,
                                     sample = HEK_samples,
                                     multisamples = "independent",
                                     plot_style = "dodge",
                                     cl = 99, colour = c("#333f50", "gray70", "#39827c", "#b5e2ff", "#ffb8c6", "#9766dc", "#ed9121"))


png(file.path(plot_dir, 'HEK_read_lengths.png'), width = 700, height = 500)
example_length_dist[["plot"]]
dev.off()


example_psite_per_region <- region_psite(reads_psite_list, annotation_df,
                                         sample = HEK_samples,
                                         plot_style = "dodge",
                                         #cl = 85,
                                         colour = c("#333f50", "gray70", "#39827c", "#b5e2ff", "#ffb8c6", "#9766dc","#ed9121"))
png(file.path(plot_dir, 'HEK_p_site_per_region.png'), width = 700, height = 500)
example_psite_per_region[["plot"]]
dev.off()

#### FILTER FOR GOOD LENGTHS ########################################################################
# no filtering for now


#### HEATMAP START/END ##############################################################################
dir_path <- file.path(plot_dir, "ribosome_start_end_heatmap")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

for (sample in names(reads_psite_list)) {
  sample_dedup_list <- list()
  sample_dedup_list[[sample]] <- c(sample)
  example_ends_heatmap <- rends_heat(reads_psite_list, 
                                     annotation_df,
                                     utr5l = 20, 
                                     cdsl = 30, 
                                     utr3l = 20, 
                                     sample = sample_dedup_list,
                                     colour = "#333f50")
  
  plot_name <- paste("plot", sample, sep = '_')
  png(file.path(plot_dir, "ribosome_start_end_heatmap",  paste(sample, 'ribosome_start_end_heatmap.png', sep = '_')),
      width = 1000, height = 500)
  print(example_ends_heatmap[[plot_name]])
  dev.off()
}



#### HEATMAP PERIODICITY ##############################################################################
# Plotting functions to check the 3 nt perioditicity
# length stratifies by read length, no length no stratification
dir_path <- file.path(plot_dir, "period_heatmap")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

for (sample in names(reads_psite_list)) {
  sample_dedup_list <- list()
  sample_dedup_list[[sample]] <- c(sample)
  example_frames_stratified <- frame_psite_length(reads_psite_list, 
                                                  annotation_df,
                                                  sample = sample_dedup_list,
                                                  region = "all",
                                                  colour = "#333f50")
  
  plot_name <- paste("plot", sample, sep = '_')
  png(file.path(plot_dir, "period_heatmap",  paste(sample, 'period_heatmap.png', sep = '_')), width = 1000, height = 500)
  print(example_frames_stratified[[plot_name]])
  dev.off()
}







####P_SITES PER FRAME ##############################################################################
dir_path <- file.path(plot_dir, "p_sites_per_frame")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

for (sample in names(reads_psite_list)) {
  sample_dedup_list <- list()
  sample_dedup_list[[sample]] <- c(sample)
  p_sites_per_frame <- frame_psite(reads_psite_list, 
                                                  annotation_df,
                                                  sample = sample_dedup_list,
                                                  region = "all",
                                                  colour = "#333f50")
  
  plot_name <- paste("plot", sample, sep = '_')
  png(file.path(plot_dir, "p_sites_per_frame",  paste(sample, 'p_sites_per_frame.png', sep = '_')), width = 1000, height = 500)
  print(p_sites_per_frame[[plot_name]])
  dev.off()
}



####METAPROFILES ##############################################################################
dir_path <- file.path(plot_dir, "metaprofiles")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

for (sample in names(reads_psite_list)) {
  sample_dedup_list <- list()
  sample_dedup_list[[sample]] <- c(sample)
  metaprofile <- metaprofile_psite(reads_psite_list, 
                                   annotation_df,
                                   sample = sample_dedup_list,
                                   utr5l = 10, 
                                   cdsl = 25, 
                                   utr3l = 10,
                                   colour = "#333f50")
  
  plot_name <- paste("plot", sample, sep = '_')
  png(file.path(plot_dir, "metaprofiles",  paste(sample, 'metaprofile.png', sep = '_')), width = 1000, height = 500)
  print(metaprofile[[plot_name]])
  dev.off()
}



####METAHEATMAP ##############################################################################
dir_path <- file.path(plot_dir, "metaheatmap")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

for (sample in names(reads_psite_list)) {
  sample_dedup_list <- list()
  sample_dedup_list[[substr(sample, 20 , nchar(sample)-6)]] <- c(sample)
  metaheatmap <- metaheatmap_psite(reads_psite_list, 
                                   annotation_df,
                                   sample = sample_dedup_list,
                                   utr5l = 10, 
                                   cdsl = 25, 
                                   utr3l = 10,
                                   colour = "#333f50")
  
  plot_name <- paste("plot", sample, sep = '_')
  png(file.path(plot_dir, "metaheatmap",  paste(sample, 'metaheatmap.png', sep = '_')), width = 1000, height = 500)
  print(metaheatmap[["plot"]])
  dev.off()
}


# there are also some codon usage plotting functions available