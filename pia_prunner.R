pia_prunner <- function(data, pia_detect_xlimit = 50, 
                        min_peak_height = 0.5, 
                        p27_pixel_cuttoff = 50,  DAPI_pixel_cutoff = 10 ) {
  
  
  # Uses the drop off in the background fluoresence you see when outside the pia with p27 and NeuN antibodies
  # calculates the point by point difference of the product of the NeuN and p27 and then determines the jump 
  # using findpeaks 
  # to manage overtrimming we exclude peaks that fall within the DAPI and p27 peak by a set value that is defined 
  # by the pixel cutoff 
  
  # p27_pixel_cuttoff: exclude pia boundaries when a p27 peak is within this number of pixels 
  # DAPI_pixel_cutoff: exclude pia boundaries when a DAPI peak is within this number of pixels 
  # min_peak_height: the size of the p27, NeuN diff product height as detected using findpeaks  
  
  require(pracma)
  require(dplyr)
  require(tidyr)
  
  
  
  
  # Prepare the data with smoothed intensity
  diff_long <- data %>%
    dplyr::filter(Row_shift_scale <= pia_detect_xlimit) %>%
    dplyr::group_by(Frame, filename, xpslit) %>%  # Group by Frame to apply smoothing within each group
    dplyr::arrange(Row_shift_scale) %>%  # Arrange by Row_shift_scale within each Frame
    dplyr::mutate(SmoothedIntensity = rollmean(RescaledIntensity, k = 5, fill = NA)) %>%
    dplyr::mutate(SmoothedIntensity = na.locf(SmoothedIntensity, na.rm = FALSE)) %>%  # Forward fill
    dplyr::mutate(SmoothedIntensity = na.locf(SmoothedIntensity, fromLast = TRUE, na.rm = FALSE)) %>%
    dplyr::mutate(diff = c(NA, diff(SmoothedIntensity))) %>%
    dplyr::mutate(diff = na.locf(diff, fromLast = TRUE, na.rm = FALSE)) %>%
    ungroup()  # Remove the grouping
  # Reshape the data to wide format
  
  all_pia_wide <- data %>%
    dplyr::select(filename, xpslit,genotype, Frame, Row_shift_scale, RescaledIntensity) %>%
    pivot_wider(id_cols =c(filename, xpslit,genotype,Row_shift_scale),names_from = Frame, values_from = RescaledIntensity, names_prefix = "scaled_intensity_") 
  
  diff_wide <- diff_long %>%
    #select(Frame, Row_shift_scale, SmoothedIntensity) %>%
    pivot_wider(id_cols = c(filename, xpslit,genotype,Row_shift_scale),
                names_from = Frame, 
                values_from = c(diff, SmoothedIntensity)) %>%
    
    dplyr::mutate(p27_NeuN_diff_product = diff_p27 * diff_NeuN) %>% 
    #dplyr::left_join(all_pia_wide, by = c("Row_shift_scale", "filename","xpslit", "genotype" )) %>%
    dplyr::mutate(pia_finder = p27_NeuN_diff_product / (SmoothedIntensity_p27 * SmoothedIntensity_NeuN +10000))
  
  
  # use find peaks to list the bumps in the p27_NeuN_diff_product
  pial_boudary_peaks <- diff_wide %>%
    dplyr::filter(Row_shift_scale <= 30) %>%
    dplyr::group_by(filename, xpslit) %>%
    dplyr::summarize(
      peak_details = list(findpeaks(p27_NeuN_diff_product, threshold = 0.3, minpeakheight = min_peak_height))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      peak_details = map(peak_details, ~ {
        # Check if the result is a matrix and has rows
        if (is.matrix(.x) && nrow(.x) > 0) {
          # Assigning names to the columns of the matrix and converting to tibble
          colnames(.x) <- c("peak_height", "peak_position", "peak_begin", "peak_end")
          as_tibble(.x)
        } else {
          # If no peaks are found or the result is not a matrix, return a tibble with NA values
          tibble(peak_height = NA_real_, peak_position = NA_integer_, peak_begin = NA_integer_, peak_end = NA_integer_)
        }
      })
    ) %>%
    tidyr::unnest(cols = peak_details)
  
  
  pial_boudary_peaks <- pial_boudary_peaks %>%
    dplyr::left_join(diff_wide, by = c("filename", "xpslit"), relationship = "many-to-many") %>%
    # Filter the rows where Row_shift_scale is within the peak range
    dplyr::filter(Row_shift_scale >= peak_begin & Row_shift_scale <= peak_end) %>%
    # Group by necessary identifiers
    dplyr::group_by(filename, xpslit, peak_height, peak_position, peak_begin, peak_end) %>%
    # Summarize the data to get the average intensities
    dplyr::summarize(
      avg_SmoothedIntensity_p27 = mean(SmoothedIntensity_p27, na.rm = TRUE),
      avg_SmoothedIntensity_NeuN = mean(SmoothedIntensity_NeuN, na.rm = TRUE),
      avg_SmoothedIntensity_DAPI = mean(SmoothedIntensity_DAPI, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
  
  
  # calculate the peak position of p27 and DAPI to exclude values 
  max_egl_tbl  <- data %>%
    # Filter to the first 100 pixels
    dplyr::filter(Row_shift_scale <= 150) %>%
    # Group by desired variables
    dplyr::group_by(Frame, filename, xpslit) %>%
    # Apply the rolling mean
    dplyr::mutate(SmoothedIntensity = zoo::rollmean(RescaledIntensity, k = 10, fill = NA)) %>%
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
  
  # Create a dataframe with all combinations of filename and xpslit
  pial_boudary_peaks_filt <- pial_boudary_peaks %>%
    #distinct(filename, xpslit) %>% 
    dplyr::left_join(max_egl_tbl_filter, by = c("filename", "xpslit")) %>%
    dplyr::filter(!(p27_exclude | DAPI_exclude)) %>%
    #dplyr::select(-p27_exclude, -DAPI_exclude ) %>%
    dplyr::filter(peak_height > min_peak_height )
  
  # First, select the peaks based on the criteria
  selected_peaks <- pial_boudary_peaks_filt %>%
    # dplyr::filter(p27_exclude | DAPI_exclude) %>%
    dplyr::group_by(filename, xpslit) %>%
    dplyr::arrange(desc(peak_height), avg_SmoothedIntensity_p27, peak_position) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Filter out peaks where avg_SmoothedIntensity_p27 is greater than 20
  filtered_peaks <- selected_peaks %>%
    dplyr::filter(avg_SmoothedIntensity_p27 <= 20)
  
  # Create a dataframe with all combinations of filename and xpslit
  all_combinations <- pial_boudary_peaks %>%
    dplyr::select(filename, xpslit) %>%
    dplyr::distinct(filename, xpslit) #%>% 
    #dplyr::left_join(max_egl_tbl, by = c("filename", "xpslit"))
  
  # Join with the filtered peaks and replace NAs with default values
  selected_peaks <- all_combinations %>%
    dplyr::left_join(filtered_peaks, by = c("filename", "xpslit")) %>%
    dplyr::mutate(peak_position = replace_na(peak_position, 0),
           peak_height = replace_na(peak_height, NA_real_),
           avg_SmoothedIntensity_p27 = replace_na(avg_SmoothedIntensity_p27, NA_real_),
           # ... same for other columns if needed
    )
  
  return(selected_peaks)
  
  
}