load_image_and_compute_rowmeans_conn2 <- function(filename, frame_names = NULL,  pixel_segment_length = 200, minimum_pixels_segment = 50) {
  # Load the image
  img <- readImage(filename)
  
  
  conn_img_list <- classify_connectivity_segments(img, small_imout = TRUE)
  conn_img <- conn_img_list$mapped_cc4
  conn_img_small <- conn_img_list$small_cc4
  
  # Extract image data
  img_data <- imageData(img)
  conn_data <- imageData(conn_img)
  conn_small_data <- imageData(conn_img_small)
  
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
      segment_data_conn <- conn_data[, segment_start:segment_end, i]
      segment_data_conn_small <- conn_small_data[, segment_start:segment_end, i]
      
      # Compute and store the mean, sd, and var for the current segment
      segment_means <- apply(segment_data, 1, mean, na.rm = TRUE)
      segment_vars <- apply(segment_data, 1, var, na.rm = TRUE)
      segment_sds <- apply(segment_data, 1, sd, na.rm = TRUE)
      
      # add the summary for the connected analysis 
      
      # ratio of background cells to other 
      segment_conn_background <- apply(segment_data_conn, 1, function(x) {
        sum(x == 0) / length(x)
      })
      
      # ratio of disconnected pixels to total 
      segment_conn_small <- apply(segment_data_conn, 1, function(x) {
        sum(x == 2) / length(x)
      })
      
      # ratio of connected pixels to total 
       segment_conn_large <- apply(segment_data_conn, 1, function(x) {
         sum(x == 1) / length(x)
       })
       
       #the number of independent cells 
       segment_conn_small_num <- apply(segment_data_conn_small, 1, function(x) {
         ((length(unique(x)) - 1) / length(x)) 
       })
       
      # Add calculated values to the main dataframe
      df$RowMean <- segment_means
      df$Variance <- segment_vars
      df$SD <- segment_sds
      df$conn_background <-  segment_conn_background
      df$conn_small <-  segment_conn_small
      df$conn_large <-  segment_conn_large
      df$conn_small_num <-  segment_conn_small_num
      
      
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