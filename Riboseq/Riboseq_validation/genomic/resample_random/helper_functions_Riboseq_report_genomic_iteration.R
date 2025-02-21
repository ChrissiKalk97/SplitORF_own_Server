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


library(glue)
library(types)
library(stringr)
library(dplyr)
library(UpSetR)
library(knitr)
library(lemon)
library(cmapR)
library(data.table)


calculate_background_threshold <- function(unique_region_type = `?`(character), path) {
    # get directory list of subdirectories of the current
    # working directory
    print(paste("region_type", unique_region_type))
    directories <- list.dirs(path = path, full.names = TRUE, recursive = TRUE)

    randomfiles <- list()
    background <- list()
    for (i in directories) {
        if (!(identical(list.files(i, pattern = paste0("*", unique_region_type,
            "_(10|[1-9])_random_intersect_counts.bed")), character(0)))) {#'|1[1-9]|20'
            # 1-20: (10|[1-9]|1[1-9]|20)
            randomfiles <- c(randomfiles, paste0(i, "/", list.files(i,
                pattern = paste0("*", unique_region_type, "_(10|[1-9])_random_intersect_counts.bed"))))
        }
    }

    randomdataframes <- lapply(randomfiles, 
                               fread, 
                               header = FALSE,
                               sep = "\t")
    
    for (randomframe in randomdataframes) {

        colnames(randomframe) <- c("chr_background", "start",
            "stop", "name", "phase", "strand", "chr_ribo", "start_ribo",
            "stop_ribo", "name_ribo", "phase_ribo", "strand_ribo",
            "nr_bp_overlap")
        
        randomframe$new_name <- paste(randomframe$name, randomframe$chr_background,
            randomframe$start, randomframe$stop, sep = "_")
        
        result <- randomframe %>%
            mutate(len = stop - start) %>%
           filter(len > 0) %>% # for the dummy regions of length 0, length 0 random regions can be sampled as well 
          group_by(new_name, name_ribo) %>%
          summarize(total_bp_overlap = sum(nr_bp_overlap, na.rm = TRUE),
                    count_combinations = n(),
                    len = first(len),
                    chr_background = first(chr_background),
                    start = first(start),
                    stop = first(stop)
                    )%>%
          # filter(total_bp_overlap > 9) %>%  # Filter based on total bp overlap
          # this is incorrect, as otherwise all 0 counts are disregarded
          group_by(new_name) %>%
          summarize(distinct_ribo_count = n_distinct(name_ribo[total_bp_overlap > 9]),
                    len = first(len),
                    chr_background = first(chr_background),
                    start = first(start),
                    stop = first(stop)) %>%
            mutate(num_reads = as.numeric(distinct_ribo_count), 
                   len = as.numeric(len),
                   relative_count = distinct_ribo_count/len) %>%
            arrange(desc(relative_count))  

        # View the result
        robust_z_scores_background <- robust_zscore(result$relative_count)
        filter_zscore_mask <- robust_z_scores_background < 3
        background_filtered <- result$relative_count[filter_zscore_mask]
        backgorund_95_percentile <- quantile(background_filtered,
            probs = c(0.95))

        background <- c(background, backgorund_95_percentile)
    }
    names(background) <- names(randomdataframes)
    return(background)
}


print_thresholds <- function(unique_region_type = `?`(character),
    background, path) {
    print(paste("region_type 2", unique_region_type))
    # get directory list of subdirectories of the current
    # working directory
    directories <- list.dirs(path = path, full.names = TRUE, recursive = TRUE)
    # obtain names of the datasets for which the thresholds
    # are calculated
    randomSetnames <- list()
    for (i in directories) {
        if (!(identical(list.files(i, pattern = paste0("*", unique_region_type,
            "_(10|[1-9])_random_intersect_counts.bed")), character(0)))) {
            randomSetnames <- c(randomSetnames, list.files(i,
                pattern = paste0("*", unique_region_type, "_(10|[1-9])_random_intersect_counts.bed")))
        }
    }
    randomnames <- stringr::str_replace(randomSetnames, pattern = "_random_intersect_counts.bed",
        replacement = "")
    # print(randomnames)

    printbackground <- data.frame(x = background)
    printbackground <- as.data.frame(t(printbackground))
    rownames(printbackground) <- randomnames
    print(randomnames)
    printbackground$sample <- stringr::str_replace(rownames(printbackground), pattern = paste0("_", unique_region_type, "_\\d+"),
                                                   replacement = "")
    printbackground <- printbackground %>% 
      group_by(sample) %>%
      summarise(threshold = mean(V1))
    rownames(printbackground) <- printbackground$sample
    print(kable(printbackground, caption = "Thresholds for each dataset"))
    return(printbackground)
}


get_intersect_files <- function(unique_region_type = `?`(character), path) {
    directories <- list.dirs(path = path, full.names = TRUE, recursive = TRUE)
    files <- list()  # Initialize an empty list
    setnames <- list()  # Initialize an empty list
    for (i in directories) {
        # Use list.files with the correct pattern to find
        # matching files
        matching_files <- list.files(i, pattern = paste0("*",
            unique_region_type, "_intersect_counts_sorted.bed"), full.names = TRUE)

        # If matching files are found, append to the lists
        if (length(matching_files) > 0) {
            files <- append(files, matching_files)
            setnames <- append(setnames, basename(matching_files))  
            # Use basename to get file names only
        }
    }


    # Read the files into a list of data frames
    dataframes <- lapply(files, fread, header = FALSE, sep = "\t")

    # Assign names to the data frames based on the setnames
    names(dataframes) <- stringr::str_replace(setnames,
                                              pattern = paste0("_", unique_region_type, "_intersect_counts_sorted.bed"),
                                              replacement = "")

    return(dataframes)
}


count_unique_regions_above_threhsold <- function(dataframes,
    threshold) {

    frames_preprocessed <- list()
    relevantregionscount <- list()
    unique_names_per_sample <- list()

    for (frame in dataframes) {
      
        colnames(frame) <- c("chr_unique", "start", "stop", "name",
            "phase", "strand", "chr_ribo", "start_ribo", "stop_ribo",
            "name_ribo", "phase_ribo", "strand_ribo", "nr_bp_overlap")


        frame$new_name <- paste(frame$name, frame$chr_unique,
            frame$start, frame$stop, sep = "_")

        processed_frame <- frame %>%
            mutate(len = stop - start) %>%
        filter(nr_bp_overlap > 0) %>%
          group_by(name, name_ribo) %>%
          summarize(total_bp_overlap = sum(nr_bp_overlap, na.rm = TRUE),
                    count_combinations = n(),
                    len = first(len),
                    chr_unique = first(chr_unique),
                    start = first(start),
                    stop = first(stop),
                    new_name = paste(new_name, collapse = "_"))%>%
          filter(total_bp_overlap > 9) %>%  # Filter based on total bp overlap
          group_by(name) %>%  # Group by name to count distinct name_ribo per name
          summarize(distinct_ribo_count = n_distinct(name_ribo),
                    len = first(len),
                    chr_unique = first(chr_unique),
                    start = first(start),
                    stop = first(stop),
                    new_name = first(new_name)) %>%
          mutate(num_reads = as.numeric(distinct_ribo_count), 
                 len = as.numeric(len),
                 relative_count = distinct_ribo_count/len) %>%
          arrange(desc(relative_count))  
          

        frames_preprocessed[[length(frames_preprocessed) +
            1]] <- processed_frame

    }


    upsetlist <- list()
    j <- 1


    for (bed in frames_preprocessed) {

        bed <- bed %>%
            filter(relative_count >= threshold[[2]][j], num_reads >
                2)
        relevantregionscount <- c(relevantregionscount, nrow(bed))
        unique_regions <- bed[, 1]
        upsetlist <- c(upsetlist, list(unique_regions))

        j <- j + 1
    }

    print(unique_regions[!(unique_regions %in% unlist(upsetlist))])

    counts_and_names_list <- list()
    counts_and_names_list$counts_above_t_relative_length <- relevantregionscount
    counts_and_names_list$unique_names <- upsetlist
    counts_and_names_list$frames_processed <- frames_preprocessed
    return(counts_and_names_list)
}


print_unique_region_above_t <- function(relevantregionscount,
    dataframes) {
    printframe <- data.frame(x = relevantregionscount)
    printframe <- as.data.frame(t(printframe))
    rownames(printframe) <- names(dataframes)
    colnames(printframe) <- paste0("Number of unique regions with relative count >= threshold")
    print(kable(printframe, caption = "Regions above the threshold",
        escape = TRUE))
}


get_top_5_unique_regions <- function(dataframes_list, threshold) {
    upsetlist <- list()
    j <- 1
    names <- names(dataframes_list)
    for (unique_region_frame in dataframes_list) {
        colnames(unique_region_frame) <-c("ID",
                                          "distinct_ribo_count",
                                          "len",
                                          "chr_unique",
                                          "start",
                                          "stop",
                                          "new_name",
                                          "num_reads",
                                          "relative_count")
        unique_region_frame %>%
            filter(relative_count >= threshold[[2]][j] & num_reads >
                2)

        if (dim(unique_region_frame)[1] > 0) {
            unique_regions <- unique_region_frame[, 1]
            upsetlist <- c(upsetlist, list(unique_regions))
            colnames(unique_region_frame) <- c("ID",
                                          "distinct_ribo_count",
                                          "len",
                                          "chr_unique",
                                          "start",
                                          "stop",
                                          "new_name",
                                          "num_reads",
                                          "relative_count")

            unique_region_frame$ID <-
            gsub('\\|', '-', unique_region_frame$ID)
            print(kable(unique_region_frame[1:5, ], caption = names(dataframes_list)[j],
                row.names = FALSE, escape = TRUE))
        } else {
            names <- names[!names %in% c(names[j])]
        }
        j <- j + 1
    }
    names(upsetlist) <- names
}


count_unique_regions_get_count <- function(dataframes, unique_names_per_sample, outdir) {
    relevantregionscount <- list()
    frame_number <- 1
    for (frame in dataframes) {
        colnames(frame) <- c("name", "num_reads", "len",
            "relative_count")
        frame$significant <- ifelse(frame$name %in% unique_names_per_sample[[frame_number]],
            1, 0)
        write.csv(frame, file.path(outdir, paste(names(dataframes)[frame_number],
            "unique_regions.csv", sep = "_")))
        above_5_counter <- 0
        c <- 1
        for (count in frame$relative_count) {
            if (frame$num_reads[[c]] >= 4 && frame$significant[[c]] == 1) {
                above_5_counter <- above_5_counter + 1
            }
            c <- c + 1
        }
        relevantregionscount <- c(relevantregionscount, above_5_counter)
        frame_number <- frame_number + 1
    }

    print <- t(as.data.frame(relevantregionscount))
    rownames(print) <- names(dataframes)
    colnames(print) <- "Count get 4"
    print(kable(print, caption = "Nr of unique regions with count get 4"))
}

Split_ORFs_validation <- function(dataframes, unique_names_per_sample, path) {
    # create empty dataframe to concatenate the dfs
    df <- data.frame(matrix(ncol = 9, nrow = 0))
    # Assign column names
    colnames(df) <- c("ID",
                    "distinct_ribo_count",
                    "len",
                    "chr_unique",
                    "start",
                    "stop",
                    "new_name",
                    "num_reads",
                    "relative_count")

    frame_number <- 1
    for (frame in dataframes) {
        colnames(frame) <- c("ID",
                            "distinct_ribo_count",
                            "len",
                            "chr_unique",
                            "start",
                            "stop",
                            "new_name",
                            "num_reads",
                            "relative_count")

        frame$significant <- ifelse(frame$ID %in% unique_names_per_sample[[frame_number]],
            1, 0)
        # print(frame)
        frame <- frame %>%
            filter(significant == 1)
        
        frame$ID <- sapply(strsplit(frame$ID , ":"), `[`, 1)

        SO_2_uniques_validated <- frame %>%
            group_by(ID) %>%
            filter(n() >= 2) %>%
            ungroup() %>%
            arrange(ID)

        write.csv(SO_2_uniques_validated[, c("ID",
                                          "distinct_ribo_count",
                                          "len",
                                          "chr_unique",
                                          "start",
                                          "stop",
                                          "new_name",
                                          "num_reads",
                                          "relative_count")], file.path(path, paste0(names(dataframes)[frame_number], '_two_regions_validated_on_transcript.csv')))
        frame_number <- frame_number + 1
        df <- rbind(df, frame)

    }
    df <- df %>%
        group_by(ID) %>%
        filter(n() >= 2) %>%
        ungroup() %>%
        arrange(ID)

}
