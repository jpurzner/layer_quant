align_folia_images2 <- function(directory_path, frame_names = NULL, reference_frame = "DAPI") {
  suppressWarnings({
    result_df_means <- load_images_from_directory(directory_path, frame_names)
  })
  
  #TODO add a default value for reference_frame 
  
  
  frame3_list <- split(result_df_means[result_df_means$Frame == reference_frame, ], result_df_means[result_df_means$Frame == reference_frame, ]$filename)
  
  # Create a dataframe to store shifts and correlation values
  alignment_df <- data.frame(reference = character(), compare = character(), shift = numeric(), correlation = numeric())
  
  # Iteratively compute shifts for each combination
  for (ref_name in names(frame3_list)) {
    for (compare_name in names(frame3_list)) {
      if (ref_name != compare_name) { # Skip comparing a curve with itself
        result <- compute_shift_and_correlation(frame3_list[[ref_name]], frame3_list[[compare_name]])
        alignment_df <- rbind(alignment_df, data.frame(reference = ref_name, compare = compare_name, shift = result$shift, correlation = result$correlation))
      }
    } # end of second for loop 
  } # end of first for loop 
  
  # Determine the file with the highest average correlation to other files
  avg_correlation <- alignment_df %>%
    dplyr::group_by(reference) %>%
    dplyr::summarize(avg_correlation = mean(correlation, na.rm = TRUE)) %>%
    dplyr::arrange(-avg_correlation)
  
  reference_filename <- as.character(avg_correlation[1, 'reference'])
  #previously we just selected the first item
  #reference_filename <- names(frame3_list)[1]
  
  # Filter alignment_df for comparisons involving the reference
  filtered_df <- alignment_df[alignment_df$reference == reference_filename, ]
  
  # Select rows with maximum correlation for each compare file
  filtered_df <- filtered_df[order(filtered_df$compare, -filtered_df$correlation), ]
  best_alignments <- filtered_df[!duplicated(filtered_df$compare), ]
  
  # Add an entry for the reference file with a shift of 0
  best_alignments <- rbind(data.frame(reference = reference_filename, 
                                      compare = reference_filename, 
                                      shift = 0, 
                                      correlation = 1), 
                           best_alignments)
  
  # Apply the shifts to the data
  aligned_data_list <- lapply(1:nrow(best_alignments), function(i) {
    align_data(best_alignments$reference[i], best_alignments$compare[i], best_alignments$shift[i], best_alignments$correlation[i], result_df_means)
  })
  # Combine into one dataframe
  aligned_data_df <- do.call(rbind, aligned_data_list)
  
return(aligned_data_df)
  }


# load all files from a directory 
load_images_from_directory <- function(start_directory, frame_names) {
  
  # Search for tif files recursively within directories named 'pieces'
  tif_files <- list.files(start_directory, pattern = "\\.tif$", recursive = TRUE,
                          include.dirs = TRUE, full.names = TRUE)
  
  # Filter out only those files within 'pieces' directories
  #tif_files <- tif_files[grep("/pieces/", tif_files)]
  
  # Apply the load_image_and_compute_rowsums function to each file
  list_of_dfs <- lapply(tif_files, function(file) load_image_and_compute_rowmeans(file, frame_names))
  
  
  # Combine the data frames, adding a 'filename' column to identify the source of each row
  combined_df <- do.call(rbind, lapply(1:length(tif_files), function(i) {
    df <- list_of_dfs[[i]]
    df$filename <- basename(tif_files[i])  # Use basename to get just the filename
    return(df)
  }))
  
  return(combined_df)
}

load_image_and_compute_rowmeans <- function(filename, frame_names = NULL) {
  # Load the image
  img <- readImage(filename)
  
  # Extract image data
  img_data <- imageData(img)
  
  # Number of frames
  num_frames <- dim(img_data)[3]
  
  if (is.null(frame_names) || length(frame_names) != num_frames) {
    frame_names <- paste0("Frame", 1:num_frames)
  }
  
  # Initialize a dataframe to store results
  df <- data.frame(Row = 1:dim(img_data)[1])
  
  # Compute row means, sd, and variance for each frame
  for (i in 1:num_frames) {
    # Compute and store the mean
    df[paste0(frame_names[i], "_Mean")] <- apply(img_data[, , i], 1, mean, na.rm = TRUE)
    
    # Compute and store the standard deviation
    df[paste0(frame_names[i], "_SD")] <- apply(img_data[, , i], 1, sd, na.rm = TRUE)
    
    # Compute and store the variance
    df[paste0(frame_names[i], "_Var")] <- apply(img_data[, , i], 1, var, na.rm = TRUE)
  }
  
  # Convert data to long format
  df_long <- df %>%
    gather(key = "Frame_Metric", value = "Value", -Row)
  
  # Separate Frame and Metric into distinct columns
  df_long <- df_long %>%
    separate(col = "Frame_Metric", into = c("Frame", "Metric"), sep = "_")
  
  return(df_long)
}

align_data <- function(ref_name, compare_name, shift_value, max_correlation, data) {
  
  # Separate data for the filename to be shifted
  compare_data <- data[data$filename == compare_name, ]
  
  aligned_data_list <- lapply(unique(compare_data$Frame), function(frame) {
    
    frame_data <- compare_data[compare_data$Frame == frame, ]
    
    # Calculate the number of values and NAs for the adjusted RowMean
    num_values <- nrow(frame_data) - abs(shift_value)
    num_nas <- nrow(frame_data) - num_values
    
    # Ensure that the sum of values and NAs matches the number of rows
    if (shift_value > 0) {
      aligned_values <- c(rep(NA, shift_value), head(frame_data$RowMean, num_values))
    } else {
      aligned_values <- c(tail(frame_data$RowMean, num_values), rep(NA, num_nas))
    }
    
    # Replace the original values with the aligned ones for the frame
    frame_data$RowMean <- aligned_values
    
    # Add the max correlation value
    frame_data$MaxCorrelation <- max_correlation
    
    return(frame_data)
    
  })
  
  # Combine the adjusted frames back together
  aligned_data <- do.call(rbind, aligned_data_list)
  
  return(aligned_data)
}


# Modified function to compute shift using ccf and also return max correlation
compute_shift_and_correlation <- function(ref_data, compare_data) {
  cc <- ccf(ref_data$RowMean, compare_data$RowMean, plot = FALSE, lag.max = length(ref_data$RowMean)/2)
  max_lag <- cc$lag[which.max(cc$acf)]
  max_correlation <- max(cc$acf)
  return(list(shift = max_lag, correlation = max_correlation))
}

