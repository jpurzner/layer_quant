egl_remove_high <- function(data) {
  
  # calculate the peak position of p27 and DAPI to exclude values 
  max_egl_tbl  <- data %>%
    # Filter to the first 100 pixels
    dplyr::filter(Row_shift_scale <= 150) %>%
    # Group by desired variables
    dplyr::group_by(Frame, filename, xpslit) %>%
    # Apply the rolling mean
    dplyr::mutate(SmoothedIntensity = zoo::rollmean(RescaledIntensityInRange, k = 10, fill = NA)) %>%
    # Identify the index of the maximum smoothed intensity
    dplyr::mutate(IndexMaxIntensity = which.max(SmoothedIntensity)) %>%
    # Use this index to get the corresponding Row_shift_scale value
    dplyr::summarize(MaxSmoothedIntensity = max(SmoothedIntensity, na.rm = TRUE),
                     RowShiftAtMax = dplyr::first(Row_shift_scale[IndexMaxIntensity])) %>%
    dplyr::ungroup() %>%
    pivot_wider(id_cols = c(filename, xpslit),
                names_from = Frame, 
                values_from = c(MaxSmoothedIntensity, RowShiftAtMax))  %>%
    dplyr::mutate(p27_exclude = ifelse(RowShiftAtMax_p27 < p27_pixel_cuttoff, TRUE, FALSE)) %>%
    dplyr::mutate(DAPI_exclude = ifelse(RowShiftAtMax_DAPI < DAPI_pixel_cutoff, TRUE, FALSE ))
  
  max_egl_tbl_filter <- max_egl_tbl %>%
    dplyr::select(filename, xpslit, p27_exclude, DAPI_exclude, RowShiftAtMax_p27, RowShiftAtMax_DAPI)
  
  
  
  
  
}