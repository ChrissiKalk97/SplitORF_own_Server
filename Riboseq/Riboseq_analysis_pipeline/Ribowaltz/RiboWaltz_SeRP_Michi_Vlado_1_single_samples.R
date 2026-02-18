# install.packages("devtools")
# library(devtools)
# install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)

library(riboWaltz)
plot_dir = '/Users/christina/Documents/own_data/SeRP/Michi_Vlado_round_1/deduplicated/Ribowaltz'
gtf_file <- '/Users/christina/Documents/own_data/Riboseq/Michi_Vlado_data_31_03_25/deduplicated/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf'
annotation_df <- create_annotation(gtfpath = gtf_file, dataSource = "RefSeqMANE0.95", organism = "Homo sapiens")


library(Biostrings) 

# Read in the FASTA file (DNA sequences example)
fasta_data <- readDNAStringSet("/Users/christina/Documents/own_data/Riboseq/Michi_Vlado_data_31_03_25/deduplicated/Ignolia_transcriptome_and_contamination.fasta")

# Retrieve the headers (sequence names)
headers <- names(fasta_data)

# Print the headers
print(headers)


####################################################
#READ in BAM FILES of deduplicated reads
####################################################


reads_list <- bamtolist(bamfolder = "/Users/christina/Documents/own_data/SeRP/Michi_Vlado_round_1/deduplicated", annotation = annotation_df)

reads_list


reads_list[["uf_muellermcnicoll_2025_04_18_IP_CHX_4_dedup "]]


psite_offset <- psite(reads_list, flanking = 6, extremity = "auto", plot = TRUE)
psite_offset

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
HEK_samples <- list('In_Puro' = c('uf_muellermcnicoll_2025_04_07_In_Puro_1_dedup',
                                       'uf_muellermcnicoll_2025_04_08_In_Puro_3_dedup',
                                       'uf_muellermcnicoll_2025_04_09_In_Puro_4_dedup'),
                      'In_CHX' = c('uf_muellermcnicoll_2025_04_10_In_CHX_1_dedup',
                                        'uf_muellermcnicoll_2025_04_11_In_CHX_2_dedup',
                                        'uf_muellermcnicoll_2025_04_12_In_CHX_4_dedup'),
                      'IP_Puro' = c('uf_muellermcnicoll_2025_04_13_IP_Puro_1_dedup',
                                    'uf_muellermcnicoll_2025_04_14_IP_Puro_3_dedup',
                                    'uf_muellermcnicoll_2025_04_15_IP_Puro_4_dedup'),
                      'IP_CHX' = c('uf_muellermcnicoll_2025_04_16_IP_CHX_1_dedup',
                                   'uf_muellermcnicoll_2025_04_17_IP_CHX_2_dedup',
                                   'uf_muellermcnicoll_2025_04_18_IP_CHX_4_dedup')
                      )


# read length distribution
example_length_dist <- rlength_distr(reads_psite_list,
                                     sample = HEK_samples,
                                     multisamples = "independent",
                                     plot_style = "dodge",
                                     cl = 99, colour = c("#333f50", "gray70", "#39827c", "#b5e2ff"))


png(file.path(plot_dir, 'HEK_read_lengths.png'), width = 700, height = 500)
example_length_dist[["plot"]]
dev.off()


example_psite_per_region <- region_psite(reads_psite_list, annotation_df,
                                         sample = HEK_samples,
                                         plot_style = "dodge",
                                         #cl = 85,
                                         colour = c("#333f50", "gray70", "#39827c"))
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