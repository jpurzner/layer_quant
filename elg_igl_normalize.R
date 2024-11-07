egl_igl_normalize <- function (data, 
                               y_name = "RowMean",
                               x_name = "Row_shift_scale",
                               igl_min_x = 250,
                               igl_max_x = 350,
                               ml_min_x = 120, 
                               ml_max_x = 160,
                               min_quant_thresh = 0.15, 
                               DAPI_max_site = "igl", 
                               smooth = TRUE) {

  
  
  
# TODO handle modify the input names as provided by y_name and x_name 
  
# smooth the data if selected  
if (smooth) {
  data <- data %>%
  dplyr::filter(Row_shift_scale <= 350) %>%
  dplyr::group_by(Frame, filename, xpslit) %>%  # Group by Frame to apply smoothing within each group
  dplyr::arrange(Row_shift_scale) %>%  # Arrange by Row_shift_scale within each Frame
  dplyr::mutate(RowMean = zoo::rollmean(RowMean, k = 10, fill = NA)) # %>%
  #dplyr:: mutate(RowMean = na.locf(RowMean, na.rm = FALSE)) %>%  # Forward fill
  #dplyr::mutate(RowMean = na.locf(RowMean, fromLast = TRUE, na.rm = FALSE))   # Backward fill
} 
  
  
  
# will normalize by predicted high and low points in the P7 cerebellum, EGL, ML and IGL.   
if (DAPI_max_site == "igl") {
  data <- data %>%
    dplyr::group_by(Frame, filename, xpslit) %>%
    dplyr::mutate(
      RescaledIntensityInRange = (RowMean - quantile(RowMean[Row_shift_scale >= 100 & Row_shift_scale <= 200], min_quant_thresh, na.rm = TRUE)) / 
          (quantile(RowMean[Row_shift_scale >= igl_min_x & Row_shift_scale <= igl_max_x], 0.95, na.rm = TRUE) - quantile(RowMean[Row_shift_scale >= ml_min_x & Row_shift_scale <= ml_max_x], min_quant_thresh, na.rm = TRUE)) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(RescaledIntensityInRange = ifelse( Frame == "DAPI",RescaledIntensityInRange / 2, RescaledIntensityInRange ))
  
  
  

  
} else if (DAPI_max_site == "egl") {
  data <- data %>%
    dplyr::group_by(Frame, filename, xpslit) %>%
    dplyr::mutate(
      RescaledIntensityInRange = case_when(
        Frame %in% c("DAPI") ~ (RowMean - quantile(RowMean[Row_shift_scale >= ml_min_x & Row_shift_scale <= ml_max_x], min_quant_thresh, na.rm = TRUE)) / 
          (quantile(RowMean[Row_shift_scale >= 50 & Row_shift_scale <= 100], 0.95, na.rm = TRUE) - quantile(RowMean[Row_shift_scale >= ml_min_x & Row_shift_scale <= ml_max_x], min_quant_thresh, na.rm = TRUE)) * 100,
        Frame %in% c("NeuN", "p27") ~ (RowMean - quantile(RowMean[Row_shift_scale >= 100 & Row_shift_scale <= 200], min_quant_thresh, na.rm = TRUE)) / 
          (quantile(RowMean[Row_shift_scale >= igl_min_x & Row_shift_scale <= igl_max_x], 0.95, na.rm = TRUE) - quantile(RowMean[Row_shift_scale >= ml_min_x & Row_shift_scale <= ml_max_x], min_quant_thresh, na.rm = TRUE)) * 100,
        TRUE ~ NA_real_  # For any unexpected cases
      )
    ) %>%
    dplyr::ungroup()
}
  
return(data)
  
}