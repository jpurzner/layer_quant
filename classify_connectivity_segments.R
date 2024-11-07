library(EBImage)
library(dplyr)

# Define the function
classify_connectivity_segments <- function(input, small_imout = FALSE, thres_imout = FALSE, pixel_cutoff = 100) {
  
  
  #TODO add optional threshold override to have specified threshold kept constant over a dataset
  
  # Check if input is a filename or an EBImage object
  if (is.character(input)) {
    img <- readImage(input, type = "tiff")
    threshold = otsu(img)
    nuc_th = EBImage::combine( mapply(function(frame, th) frame > th, getFrames(img), threshold, SIMPLIFY=FALSE) )
    #TODO add a another else if to manage 
  } else if (inherits(input, "Image")) {
    nuc_th <- input
  } else {
    stop("Input must be a filename or an EBImage object.")
  }

  
  # if no threshold is provided then it will be calculated from each image
 
  cc4 <- bwlabel(nuc_th)
  
  # Initialize an empty data frame to store the results
  pixel_counts_df_all_frames <- data.frame()
  
  # Number of frames in cc4
  num_frames <- dim(cc4)[3]
  
  # Iterate through each frame
  for (frame in 1:num_frames) {
    # Get pixel counts for the current frame
    pixel_counts <- table(cc4[,,frame])
    
    # Convert to a data frame
    pixel_counts_df <- as.data.frame(pixel_counts)
    names(pixel_counts_df) <- c("Segment", "Count")
    
    # Process the data frame
    pixel_counts_df <- pixel_counts_df %>%
      dplyr::mutate(Segment = as.numeric(as.character(Segment)),  # Convert Segment from factor to numeric
             Frame = frame) %>%  # Add a column for the frame number
      dplyr::arrange(desc(Count)) %>%
      dplyr::mutate(Rank = row_number()) %>%
      dplyr::mutate(mapsegment = case_when(
        Segment == 0 ~ 0,
        Segment > 0 & Count > pixel_cutoff ~ 1,
        Segment > 0 & Count <= pixel_cutoff ~ 2,
        TRUE ~ NA_integer_  # This line handles any unexpected cases
      )) %>%
        dplyr::mutate(smallonly = ifelse(mapsegment ==2, Segment, 0 )) 
     
    
    # Combine the results for each frame
    pixel_counts_df_all_frames <- rbind(pixel_counts_df_all_frames, pixel_counts_df)
  }
  
  mapped_cc4 <- cc4
  small_mapped_cc4 <- cc4 
  
  # Iterate through each frame
  for (frame in 1:num_frames) {
    # Get the mapping for the current frame
    mapping <- pixel_counts_df_all_frames %>%
      dplyr::filter(Frame == frame) %>%
      dplyr::select(Segment, mapsegment)
    
    mapping_small <- pixel_counts_df_all_frames %>%
      dplyr::filter(Frame == frame) %>%
      dplyr::select(Segment, smallonly)
    
    # Convert to a named vector for easy mapping
    mapping_vector <- setNames(mapping$mapsegment, mapping$Segment)
    small_mapping_vector <- setNames(mapping_small$smallonly, mapping_small$Segment)
    
    # Map the values
    mapped_cc4[,,frame] <- mapping_vector[as.character(cc4[,,frame])]
    small_mapped_cc4[,,frame] <- small_mapping_vector[as.character(cc4[,,frame])]
  }
  
  
  # Add to the return logic
  result_list <- list(mapped_cc4 = mapped_cc4)
  
  if (small_imout) {
    result_list$small_cc4 <- small_mapped_cc4
  }
  
  if (thres_imout) {
    result_list$nuc_th <- nuc_th
  }
  
  # Return results based on the options provided
  if (length(result_list) == 1) {
    # If there is only one element in the list, return it directly
    return(result_list$mapped_cc4)
  } else {
    # Otherwise, return the list of results
    return(result_list)
  }
}
  

