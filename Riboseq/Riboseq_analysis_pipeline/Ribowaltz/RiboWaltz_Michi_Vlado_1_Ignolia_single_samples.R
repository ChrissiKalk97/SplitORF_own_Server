# install.packages("devtools")
# library(devtools)
# install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)

library(riboWaltz)
library(Biostrings)

plot_dir = '/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/Ribowaltz'
gtf_file <- '/projects/splitorfs/work/Riboseq/data/contamination/Ignolia_paper/mRNA/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf'
annotation_df <- create_annotation(gtfpath = gtf_file, dataSource = "RefSeqMANE0.95", organism = "Homo sapiens")

 

# Read in the FASTA file (DNA sequences example)
fasta_data <- readDNAStringSet("/projects/splitorfs/work/reference_files/own_data_refs/Riboseq/Ignolia/Ignolia_transcriptome_and_contamination.fasta")

# Retrieve the headers (sequence names)
headers <- names(fasta_data)



####################################################
#READ in BAM FILES of deduplicated reads
####################################################


reads_list <- bamtolist(bamfolder = "/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_concat_transcriptome_Ignolia/filtered/q10/dedup", annotation = annotation_df)



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

#### FILTER FOR GOOD LENGTHS ########################################################################
reads_psite_list <- length_filter(data = reads_psite_list,
                               length_filter_mode = "custom",
                               length_range = 33:40)


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