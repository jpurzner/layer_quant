plot_mcp_segments2 <- function(filename_example, sign, channel,  means_all, data_pruned_shift_rescaled_split_all) {
  
  # Filter the data for the given filename and xpslit
  test_image <- data_pruned_shift_rescaled_split_all %>%
    dplyr::filter(filename == filename_example) %>%
    dplyr::filter(xpslit == sign) %>%
    dplyr::filter(Frame == channel)
  
  # Extract relevant values from the means_all dataframe
  extracted_values <- means_all %>%
    dplyr::filter(filename == filename_example & xpslit == sign) 
  
  # Extract change points and slopes from extracted_values
  change_points <- as.numeric(extracted_values[, grepl("cp_", names(extracted_values))])
  slopes <- as.numeric(extracted_values[, grepl("Row_shift_scale_", names(extracted_values))])
  intercept <- as.numeric(extracted_values$int_1)
  
  # Handling NA in slopes:
  na_index <- is.na(slopes)
  slopes <- slopes[!na_index]
  # If it affects change_points, adapt them as well:
  if(length(change_points) > length(slopes)){
    change_points <- change_points[-length(change_points)] 
  }
  
  intercepts <- numeric(length(slopes))
  y_cp <- numeric(length(change_points))
  
  intercepts[1] <- intercept
  for (i in 1:(length(slopes) - 1)) {
    y_cp[i] <- slopes[i] * change_points[i] + intercepts[i]
    intercepts[i + 1] <- y_cp[i] - slopes[i + 1] * change_points[i]
  }
  
  # Setup the ggplot with the data
  p <- ggplot(test_image, aes(x = Row_shift_scale, y = RescaledIntensity)) +
    geom_line() +
    geom_vline(xintercept = as.numeric(change_points), linetype="dashed", color = "red") +
    labs(title = "Changepoint Analysis",
         x = "Observation",
         y = "Rescaled Intensity") +
    theme_minimal()
  
  # Add segments to the plot
  # Add segments to the plot
  for (i in 1:length(slopes)) {
    if (i == 1) {
      xstart = min(test_image$Row_shift_scale)
      xend = change_points[i]
    } else if (i == length(slopes)) {
      xstart = change_points[i-1]
      xend = max(test_image$Row_shift_scale)
    } else {
      xstart = change_points[i-1]
      xend = change_points[i]
    }
    # Calculate the corresponding y values for the start and end of the segment
    ystart = slopes[i] * xstart + intercepts[i]
    yend = slopes[i] * xend + intercepts[i]
    
    # Adding the segment using geom_line
    p <- p + geom_line(data = data.frame(x = c(xstart, xend), y = c(ystart, yend)),
                       aes(x = x, y = y), color = "blue")
  }
  
  return(p)
}




