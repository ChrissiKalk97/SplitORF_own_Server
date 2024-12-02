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


library(glue)
library(types)
library(stringr)
library(dplyr)
library(UpSetR)
library(knitr)
library(lemon)
library(cmapR)


calculate_background_threshold <- function(unique_region_type = `?`(character), path) {
    # get directory list of subdirectories of the current
    # working directory
    directories <- list.dirs(path = path, full.names = TRUE, recursive = TRUE)

    randomfiles <- list()
    background <- list()
    for (i in directories) {
        if (!(identical(list.files(i, pattern = paste0("*", unique_region_type,
            "_random_intersect_counts.bed")), character(0)))) {
            randomfiles <- c(randomfiles, paste0(i, "/", list.files(i,
                pattern = paste0("*", unique_region_type, "_random_intersect_counts.bed"))))
        }
    }
    print(randomfiles)

    randomdataframes <- lapply(randomfiles, read.csv, header = FALSE,
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
            group_by(new_name) %>%
            summarize(num_reads = ifelse(sum(nr_bp_overlap >
                0) > 0, sum(nr_bp_overlap > 0), 0), len = first(len)) %>%
            filter(len > 0) %>%
            mutate(num_reads = as.numeric(num_reads), len = as.numeric(len),
                relative_count = num_reads/len) %>%
            arrange(desc(relative_count))  

        # View the result
        print(result)
        robust_z_scores_background <- robust_zscore(result$relative_count)
        filter_zscore_mask <- robust_z_scores_background < 3
        background_filtered <- result$relative_count[filter_zscore_mask]
        backgorund_95_percentile <- quantile(background_filtered,
            probs = c(0.95))
        print(backgorund_95_percentile)

        background <- c(background, backgorund_95_percentile)
    }
    return(background)
}


print_thresholds <- function(unique_region_type = `?`(character),
    background, path) {
    # get directory list of subdirectories of the current
    # working directory
    directories <- list.dirs(path = path, full.names = TRUE, recursive = TRUE)
    # obtain names of the datasets for which the thresholds
    # are calculated
    randomSetnames <- list()
    for (i in directories) {
        if (!(identical(list.files(i, pattern = paste0("*", unique_region_type,
            "_random_intersect_counts.bed")), character(0)))) {
            randomSetnames <- c(randomSetnames, list.files(i,
                pattern = paste0("*", unique_region_type, "_random_intersect_counts.bed")))
        }
    }
    randomnames <- stringr::str_replace(randomSetnames, pattern = "_random_intersect_counts_relative_sorted.bed",
        replacement = "")

    printbackground <- data.frame(x = background)
    printbackground <- as.data.frame(t(printbackground))
    rownames(printbackground) <- randomnames
    colnames(printbackground) <- "Threshold"
    print(kable(printbackground, caption = "Thresholds for each dataset"))
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
    dataframes <- lapply(files, read.csv, header = FALSE, sep = "\t")

    # Assign names to the data frames based on the setnames
    names(dataframes) <- stringr::str_replace(setnames,
                                              pattern = "_intersect_counts_sorted.bed",
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
            filter(len > 0) %>%
            group_by(new_name) %>%
            summarize(num_reads = ifelse(sum(nr_bp_overlap >
                0) > 0, sum(nr_bp_overlap > 0), 0), len = first(len)) %>%
            mutate(num_reads = as.numeric(num_reads), len = as.numeric(len),
                relative_count = num_reads/len) %>%
            arrange(desc(relative_count))

        frames_preprocessed[[length(frames_preprocessed) +
            1]] <- processed_frame

    }


    upsetlist <- list()
    j <- 1
    i <- 1

    for (bed in frames_preprocessed) {

        bed <- bed %>%
            filter(relative_count >= threshold[[j]], num_reads >
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
        colnames(unique_region_frame) <- c("new_name", "num_reads",
            "len", "relative_count")
        unique_region_frame %>%
            filter(relative_count >= threshold[[j]] & num_reads >
                2)

        if (dim(unique_region_frame)[1] > 0) {
            unique_regions <- unique_region_frame[, 1]
            upsetlist <- c(upsetlist, list(unique_regions))
            colnames(unique_region_frame) <- c("ID", "num_reads",
                "len", "relative_count")

            # unique_region_frame$ID <-
            # paste0(str_split_fixed(unique_region_frame$ID,
            # '\\|', 2)[,1], '-',
            # str_split_fixed(unique_region_frame$ID,
            # '\\|', 2)[,2])
            print(kable(unique_region_frame[1:5, ], caption = names(dataframes_list)[j],
                row.names = FALSE, escape = TRUE))
        } else {
            names <- names[!names %in% c(names[j])]
        }
        j <- j + 1
    }
    names(upsetlist) <- names
    return(upsetlist)
}


count_unique_regions_get_count <- function(dataframes, unique_names_per_sample, outdir) {
    relevantregionscount <- list()
    frame_number <- 1
    for (frame in dataframes) {
        colnames(frame) <- c("new_name", "num_reads", "len",
            "relative_count")
        frame$significant <- ifelse(frame$new_name %in% unique_names_per_sample[[frame_number]],
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

Split_ORFs_validation <- function(dataframes, unique_names_per_sample) {
    # create empty dataframe to concatenate the dfs
    df <- data.frame(matrix(ncol = 6, nrow = 0))
    # Assign column names
    colnames(df) <- c("new_name", "num_reads", "len",
                      "relative_count")

    frame_number <- 1
    for (frame in dataframes) {
        colnames(frame) <- c("new_name", "num_reads", "len",
                             "relative_count")
        frame$significant <- ifelse(frame$new_name %in% unique_names_per_sample[[frame_number]],
            1, 0)
        # print(frame)
        frame <- frame %>%
            filter(significant == 1)
        
        frame$ID <- sapply(strsplit(frame$new_name, ":"), `[`, 1)

        SO_2_uniques_validated <- frame %>%
            group_by(ID) %>%
            filter(n() >= 2) %>%
            ungroup() %>%
            arrange(ID)

        print(SO_2_uniques_validated[, c("new_name", "num_reads", "len",
                                         "relative_count")], max = nrow(SO_2_uniques_validated))
        frame_number <- frame_number + 1
        df <- rbind(df, frame)

    }
    df <- df %>%
        group_by(ID) %>%
        filter(n() >= 2) %>%
        ungroup() %>%
        arrange(ID)


    # print(df[, c("new_name", "num_reads", "len",
      #           "relative_count")],
       # max = nrow(df))
}
