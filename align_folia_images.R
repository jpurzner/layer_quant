align_folia_images <- function(directory_path, frame_names = NULL, reference_frame = "DAPI", pixel_segment_length = 200) {
  
  # Check if EBImage is installed; if not, install it
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!requireNamespace("EBImage", quietly = TRUE))
    BiocManager::install("EBImage")
  
  # Load the necessary library
  library(EBImage)
  
  suppressWarnings({
    result_df_means <- load_images_from_directory(directory_path, frame_names, pixel_segment_length)
  })
  
  # Create a list of data frames, one for each unique filename in the reference frame
  reference_list <- split(result_df_means[result_df_means$Frame == reference_frame, ], 
                          result_df_means[result_df_means$Frame == reference_frame, ]$filename)
  
  alignment_df <- data.frame(reference = character(), compare = character(), shift = numeric(), correlation = numeric())
  
  # Iteratively compute shifts for each combination
  for (ref_name in names(reference_list)) {
    for (compare_name in names(reference_list)) {
      if (ref_name != compare_name) { # Skip comparing a curve with itself
        result <- compute_shift_and_correlation(reference_list[[ref_name]], reference_list[[compare_name]])
        alignment_df <- rbind(alignment_df, data.frame(reference = ref_name, compare = compare_name, shift = result$shift, correlation = result$correlation))
      }
    }
  }
  
  # Determine the file with the highest average correlation to other files
  avg_correlation <- alignment_df %>%
    dplyr::group_by(reference) %>%
    dplyr::summarize(avg_correlation = mean(correlation, na.rm = TRUE)) %>%
    dplyr::arrange(-avg_correlation)
  
  # The top reference is the one with the highest average correlation
  top_reference <- avg_correlation[1, ]$reference
  
  
  
  # Filter alignment_df for those rows where the reference is the top_reference
  alignment_to_apply <- alignment_df[alignment_df$reference == top_reference, ]
  
  # Apply the shifts to the data and align
  aligned_data_list <- lapply(1:nrow(alignment_to_apply), function(i) {
    align_data(alignment_to_apply$reference[i], alignment_to_apply$compare[i], alignment_to_apply$shift[i], alignment_to_apply$correlation[i], result_df_means)
  })
  
  # Combine into one dataframe
  aligned_data_df <- do.call(rbind, aligned_data_list)
  
  return(aligned_data_df)
}

  
  
  # load all files from a directory 
load_images_from_directory <- function(start_directory, frame_names, pixel_segment_length) {
  
  # Search for tif files recursively within the directory
  tif_files <- list.files(start_directory, pattern = "\\.tif$", recursive = TRUE,
                          include.dirs = TRUE, full.names = TRUE)
  
  # Apply the load_image_and_compute_rowmeans function to each file
  list_of_dfs <- lapply(tif_files, function(file) load_image_and_compute_rowmeans(file, frame_names, pixel_segment_length))
  
  # Combine the data frames, adding a 'filename' column to identify the source of each row
  combined_df <- do.call(rbind, list_of_dfs)
  return(combined_df)
}

  
  load_image_and_compute_rowmeans <- function(filename, frame_names = NULL, pixel_segment_length, minimum_pixels_segment = 50) {
    # Load the image
    img <- readImage(filename)
    
    # Extract image data
    img_data <- imageData(img)
    
    # Number of frames
    num_frames <- dim(img_data)[3]
    
    if (is.null(frame_names) || length(frame_names) != num_frames) {
      frame_names <- paste0("Frame", 1:num_frames)
    }
    
    # Calculate the number of segments based on the provided pixel segment length
    num_cols <- dim(img_data)[2]
    full_segments <- floor(num_cols / pixel_segment_length) # The number of full segments
    remainder_pixels <- num_cols %% pixel_segment_length # The number of pixels in the partial segment
    
    # Check if the remainder meets the minimum pixel requirement
    if (remainder_pixels >= minimum_pixels_segment) {
      num_segments <- full_segments + 1
    } else {
      num_segments <- full_segments
    }
    
    # Number of rows (preserved across segments)
    num_rows <- dim(img_data)[1]
    
    # Initialize a dataframe to store results
    df <- data.frame(Row = 1:num_rows) # Segment column removed
    
    # Empty dataframe to store long-format data
    df_long <- data.frame()
    
    # Compute row means, sd, and variance for each frame
    for (i in 1:num_frames) {
      frame_name <- ifelse(is.null(frame_names), paste0("Frame", i), frame_names[i])
      
      # Split the data into segments and compute the metrics for each segment
      for (j in 1:num_segments) {
        segment_start <- (j - 1) * pixel_segment_length + 1
        segment_end <- ifelse(j == num_segments & remainder_pixels >= minimum_pixels_segment, 
                              num_cols, 
                              j * pixel_segment_length)
        
        # Select the segment of the image data
        segment_data <- img_data[, segment_start:segment_end, i]
        
        # Compute and store the mean, sd, and var for the current segment
        segment_means <- apply(segment_data, 1, mean, na.rm = TRUE)
        segment_vars <- apply(segment_data, 1, var, na.rm = TRUE)
        segment_sds <- apply(segment_data, 1, sd, na.rm = TRUE)
        
        # Add calculated values to the main dataframe
        df$RowMean <- segment_means
        df$Variance <- segment_vars
        df$SD <- segment_sds
        
        # Add a column with the frame name
        df$Frame <- frame_name
        
        # Modify the filename to include segment information
        filename_base <- gsub("\\.tif$", "", basename(filename))  # remove .tif from the filename
        modified_filename <- paste0(filename_base, "_segment_", j) # append segment number
        
        # Adding modified filename to the dataframe
        df$filename <- modified_filename
        
        # Binding the frames in long format
        df_long <- rbind(df_long, df)
      }
    }
    
    return(df_long)
  }
  
  
  
  
  
  align_data <- function(ref_name, compare_name, shift_value, max_correlation, data) {
    # Filter data for the specific filename to be shifted
    compare_data <- data[data$filename == compare_name, ]
    
    # Extract columns that contain metric data
    metric_cols <- c("RowMean", "Variance", "SD")
    
    # Initialize an empty data frame for aligned data, including the 'Frame' column
    aligned_data <- data.frame(Row = compare_data$Row, Frame = compare_data$Frame)  # Include Frame here
    
    # Process each metric column
    for (metric_col in metric_cols) {
      metric_data <- compare_data[, metric_col]
      
      # Calculate the number of values and NAs for the adjusted Value
      num_values <- nrow(compare_data) - abs(shift_value)
      num_nas <- nrow(compare_data) - num_values
      
      # Ensure that the sum of values and NAs matches the number of rows
      if (shift_value > 0) {
        aligned_values <- c(rep(NA, shift_value), head(metric_data, num_values))
      } else {
        aligned_values <- c(tail(metric_data, num_values), rep(NA, num_nas))
      }
      
      # Add the aligned metric to the aligned_data
      aligned_data[[metric_col]] <- aligned_values
    }
    
    # Assign max correlation
    aligned_data$MaxCorrelation <- max_correlation
    
    # Maintain filename for traceability
    aligned_data$filename <- compare_name
    
    # The line for including the segment in the output is removed, as segments are no longer used.
    
    return(aligned_data)
  }
  
  
  
  
  compute_shift_and_correlation <- function(ref_data, compare_data) {
    # Initialize max_lag and max_correlation to 0
    max_lag <- 0
    max_correlation <- 0
    
    # Check if either dataset is empty or contains only NA values
    if (nrow(ref_data) == 0 || nrow(compare_data) == 0 ||
        all(is.na(ref_data$RowMean)) || all(is.na(compare_data$RowMean))) {
      warning("Insufficient data - data is empty or contains only NAs.")
      return(list(shift = max_lag, correlation = max_correlation)) # returns 0s in case of problematic data
    }
    
    # Check for sufficient number of observations
    min_observations <- 10 # set a threshold for the minimum acceptable number of observations
    if (nrow(ref_data) < min_observations || nrow(compare_data) < min_observations) {
      warning(paste("Not enough observations - need at least", min_observations, "observations."))
      return(list(shift = max_lag, correlation = max_correlation)) # returns 0s in case of insufficient observations
    }
    
    # Try to compute ccf and catch any potential errors
    tryCatch({
      cc <- ccf(ref_data$RowMean, compare_data$RowMean, plot = FALSE, lag.max = length(ref_data$RowMean)/2)
      max_lag <- cc$lag[which.max(cc$acf)]
      max_correlation <- max(cc$acf)
    }, error = function(e) {
      warning(paste("Error computing cross-correlation:", e$message))
      # returns 0s in case of an error during computation
    })
    
    return(list(shift = max_lag, correlation = max_correlation))
  }
  
  