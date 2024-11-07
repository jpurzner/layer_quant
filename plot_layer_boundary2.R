plot_layer_boundary2 <- function(data, extra_anno = NULL, 
                                 sel_frame = "p27", 
                                 window_size = 5, 
                                 pixel_width = NULL, 
                                 x_col_name = "Row_shift_scale_micron",
                                 x_type = "micron", 
                                 x_min = 0,
                                 x_max = 150,
                                 y_min = 0, 
                                 y_max = 150,
                                 y_col_name = "RescaledIntensity",
                                 comparison_name = "genotype", 
                                 comparison_order = c("WT", "Ezh2_cKO"), 
                                 remove_underscore = FALSE, 
                                 metadata = extra_anno, 
                                 biological_rep = NULL, 
                                 label_reps = FALSE, 
                                 stat = "ttest") {

  
# data must contain filename and xpslit values from the upstream processing.   
  
# label_reps: when TRUE the biological replicates are labeled using ggrepel   
  
# Examples 
# plot_layer_boundary2(data_pia_aligned_micron, x_col_name = "Row_shift_scale_pia_micron", y_col_name = "RescaledIntensityInRange", comparison_order = c("WT", "Ezh2_cKO"))  
  
# TODO fix the data trimming for out of bounds 

  
  
# load libraries 
require(patchwork)
require(lmerTest)
require(RColorBrewer)
require(ggplot2)
require(dplyr)
require(ggsignif)

#rename the data columns to make data manipulation easier 
names(data)[names(data) == x_col_name] <- "row_index"
names(data)[names(data) == comparison_name] <- "comparison"  
names(data)[names(data) == y_col_name] <- "intensity"    
  
# set up comparison column 
if(!is.factor(data$comparison)) {
  if (!is.null(comparison_order)) {
    data$comparison = factor(data$comparison, levels = comparison_order)
  } else {
    data$comparison <- factor(data$comparison)   
  } 
}

# convert pixels to microns if value provided
# for the original data set the pixel width was 
if (!is.null(pixel_width)) {
  data$row_index = data$row_index * pixel_width
}


comparison_levels <- levels(data$comparison)
  

# check the provided data is consistend with options 
 

# set variables for colour   
genotype_colors <- brewer.pal(2, "Set1")
color_mapping <- setNames(genotype_colors, c("WT", "Ezh2_cKO"))
dashed_color <- brewer.pal(4, "Set1")[4]

# calculate the average for each Row_shift_scale within each facet
averages <- data %>%
  dplyr::ungroup() %>%
  dplyr::filter(Frame == sel_frame) %>%
  dplyr::group_by(row_index, comparison) %>%
  dplyr::summarize(average_intensity = mean(intensity, na.rm = TRUE)) 

averages_wide <- averages %>%
  pivot_wider(names_from = comparison, values_from = average_intensity) %>%
  # removed redundant data on the lower and upper limit 
  dplyr::filter(row_index > x_min & row_index < x_max) %>%
  # Calculate the difference between genotypes
  dplyr::mutate(Difference = !!sym(comparison_levels[2]) - !!sym(comparison_levels[1])) 


# calculate the max and min change values
# TODO apply limits to where the max and min can be found 
max_diff_index_pos <- which.max(averages_wide$Difference)
min_diff_index_pos <- which.min(averages_wide$Difference)

# Retrieve the corresponding Row_shift_scale_micron
max_diff_row_index <- averages_wide$row_index[max_diff_index_pos]
min_diff_row_index <- averages_wide$row_index[min_diff_index_pos]

# Set up the windows 
# Need to allow for non integer values as user may provide micron values 
# Calculate the start and end indices
half_n <- window_size %/% 2
start_index <- max(1, max_diff_index_pos - half_n)
end_index <- min(length(averages_wide$row_index), max_diff_index_pos + half_n)
max_diff_row_index_window <- averages_wide$row_index[start_index:end_index]
max_diff_row_index_window <- data.frame(TargetRow = max_diff_row_index, row_index = max_diff_row_index_window)

half_n <- window_size %/% 2
start_index <- max(1, min_diff_index_pos - half_n)
end_index <- min(length(averages_wide$row_index), min_diff_index_pos + half_n)
min_diff_row_index_window <- averages_wide$row_index[start_index:end_index]
min_diff_row_index_window <- data.frame(TargetRow = min_diff_row_index, row_index = min_diff_row_index_window)

expanded_ranges <-rbind(max_diff_row_index_window, min_diff_row_index_window)


# generate the line plot 
frame_line <- data %>%
  filter(Frame == sel_frame) %>%
  ggplot() + 
  geom_line(aes(x = row_index, y = intensity, group  = paste0(filename, xpslit)), 
            colour = "black", alpha = 0.05) +
  geom_line(data = averages, aes(x = row_index, y = average_intensity, group = comparison, color = comparison), 
            size = 1) +
  geom_vline(xintercept = max_diff_row_index, linetype = "dashed", color = dashed_color) +
  geom_vline(xintercept = min_diff_row_index, linetype = "dashed", color = dashed_color) +
  scale_color_manual(values = color_mapping) +
  labs(title = "",
       x = "",
       y = paste("Normalized", sel_frame, "\nmean intensity")) +  # Dynamic y-axis label
  xlim(x_min, x_max) +
  ylim(y_min, y_max) +
  facet_grid(comparison ~ Frame) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = -1, unit = "cm")  # Negative bottom margin to bring plots closer
  )


# Plot the difference
avg_difference_plot <- ggplot(averages_wide, aes(x = row_index, y = Difference)) +
  geom_line() +
  geom_vline(xintercept = max_diff_row_index, linetype = "dashed", color = dashed_color) +
  geom_vline(xintercept = min_diff_row_index, linetype = "dashed", color = dashed_color) +
  scale_x_continuous(position = "top") +  # Move x-axis to top
  labs(title = "",
       x = "Micrometers from pial boundary",
       y = "Intensity Difference\nEzh2 cKO - WT") +
  xlim(0, 150) +
  theme_minimal() + 
  theme(
    plot.margin = margin(t = -1, unit = "cm"),  # Adjust top margin
    legend.position = "none"
  )

target_rows <- c( max_diff_index_pos, min_diff_index_pos)

# Merge expanded_ranges with the filtered data
filtered_expanded_data <- data %>%
  filter(Frame == sel_frame) %>%
  inner_join(expanded_ranges, by = "row_index")

averaged_data <- filtered_expanded_data %>%
  dplyr::group_by(filename, xpslit, comparison, TargetRow) %>%
  dplyr::summarize(AverageIntensity = mean(intensity, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Merge with the original data to get the associated metadata
merged_data <- dplyr::left_join(averaged_data, 
                                data %>% 
                                  dplyr::select(filename, xpslit, comparison, row_index, dataset), 
                                by = c("filename", "xpslit", "comparison", "TargetRow" = "row_index"))

# Ensure unique rows by selecting the first occurrence of each group
final_data <- merged_data %>%
  dplyr::group_by(filename, xpslit, comparison, TargetRow) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

if (!is.null(metadata)) {
  # check that all additional columns are factors   
  final_data <- dplyr::left_join(final_data, metadata, by = c("filename", "xpslit"))
  
}

# biological replicates are provided they will be shown as a dots within the violin plot 
if (!is.null(biological_rep)) {
  #check if the biological replicates are within the data 
  
  # convert names so easier to manage 
  names(final_data)[names(final_data) == biological_rep] <- "biological_rep"
  
  # average across the biological replicates 
  bio_rep_average <- final_data %>%
    dplyr::group_by(biological_rep, comparison, TargetRow) %>%
    dplyr::summarize(AverageIntensity = mean(AverageIntensity, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  ##CALCULATE THE PVALUE 
  # split of the max 
  max_diff_final_data_ano <- final_data %>%
    dplyr::filter(TargetRow == max_diff_row_index)  
  max_diff_bio_rep <- bio_rep_average %>%
    dplyr::filter(TargetRow == max_diff_row_index)
  # calculate the p-values 
  if(stat == "lmer") {
    model <- lmer(AverageIntensity ~ comparison +  ( 1 | biological_rep ), data = max_diff_final_data_ano)
    model_summary <- summary(model)
    p_values_max <- model_summary$coefficients[, "Pr(>|t|)"] 
    p_values_max <- p_values_max[2] 
  } else if (stat == "ttest") {
    model_summary <- t.test(AverageIntensity ~ comparison, data = max_diff_bio_rep)  
    p_values_max <- model_summary$p.value
  }
  
  # split off the min
  min_diff_final_data_ano <- final_data %>%
    dplyr::filter(TargetRow == min_diff_row_index)    
  min_diff_bio_rep <- bio_rep_average %>%
    dplyr::filter(TargetRow == min_diff_row_index)
  if(stat == "lmer") {
    model <- lmer(AverageIntensity ~ comparison +  ( 1 | biological_rep ), data = min_diff_final_data_ano)
    model_summary <- summary(model)
    p_values_min <- model_summary$coefficients[, "Pr(>|t|)"]
    p_values_min <- p_values_min[2] 
  } else if (stat == "ttest") {
    model_summary <- t.test(AverageIntensity ~ comparison, data = min_diff_bio_rep)  
    p_values_min <- model_summary$p.value
  }
  
  #micron_min <- unique(final_data[final_data$TargetRow == min_diff_row_index, "Row_shift_scale_micron"])
  #micron_max <- unique(final_data[final_data$TargetRow == max_diff_row_index, "Row_shift_scale_micron"])
  
  # Create custom titles
  title_min <- paste("Micrometers from\npial boundary:", round(min_diff_row_index, 2))
  title_max <- paste("Micrometers from\npial boundary:", round(max_diff_row_index, 2))
  centered_title_theme <- theme(plot.title = element_text(hjust = 0.5, size = 12))  # Adjust the size as needed
  
  #TODO need to handle the case when the max and min are backwards, should be first and second as it appears in space rather then max or min 
  
  if (label_reps) {
    # when label_reps is TRUE ggrepel is used to label the replicates
    
      pos <- position_jitter(width = 0.5, seed = 1)
      
      plot_min <- ggplot(subset(final_data, TargetRow == min_diff_row_index), aes(x = comparison, y = AverageIntensity, fill = comparison)) +
      geom_violin(trim = FALSE) +
      geom_signif(comparisons = list(c("Ezh2_cKO", "WT")), annotations = paste("p =", formatC(p_values_min, format = "e", digits = 2)), y_position = max(final_data$AverageIntensity) * 0.9) +
      geom_point(data = subset(bio_rep_average, TargetRow == min_diff_row_index),  
                  aes(x = comparison, y = AverageIntensity, group = comparison), position = pos, width = 0.2, size = 1, colour = "black") +
      geom_text_repel(data = subset(bio_rep_average, TargetRow == min_diff_row_index),
                       aes(x = comparison, y = AverageIntensity, group = comparison, label = biological_rep),
                       position =pos,
                         #box.padding = 0.35, 
                         #point.padding = 0.5,
                         segment.color = 'black') +
      labs(x = "Genotype", 
           y = paste("Normalized", sel_frame, "\nmean intensity")) +  # Dynamic y-axis label
      ylim(-10, 150) +
      scale_fill_brewer(palette = "Set1") +
      labs(title = title_min) +
      theme_minimal() + 
      centered_title_theme +
      theme(
        legend.position = "none")
    
      plot_max <- ggplot(subset(final_data, TargetRow == max_diff_row_index), aes(x = comparison, y = AverageIntensity, fill = comparison)) +
      geom_violin(trim = FALSE) +
      geom_signif(comparisons = list(c("Ezh2_cKO", "WT")), annotations = paste("p =", formatC(p_values_max, format = "e", digits = 2)), y_position = max(final_data$AverageIntensity) * 0.9) +
      geom_point(data = subset(bio_rep_average, TargetRow == max_diff_row_index),  
                 aes(x = comparison, y = AverageIntensity, group = comparison), 
                 position = pos, 
                 width = 0.2, 
                 size = 1, 
                 colour = "black") +
      geom_text_repel(data = subset(bio_rep_average, TargetRow == max_diff_row_index),
                       aes(x = comparison, y = AverageIntensity, group = comparison, label = biological_rep),
                        position =pos,
                         #box.padding = 0.35, 
                         #point.padding = 0.5,
                         segment.color = 'black') +
      labs(x = "Genotype", 
           y = paste("Normalized", sel_frame, "\nmean intensity")) +  # Dynamic y-axis label
      scale_fill_brewer(palette = "Set1") +
      ylim(-10, 150) +
      labs(title = title_max) +
      theme_minimal() +
      centered_title_theme +
      theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    
    
  } else {
    # when label_reps is FALSE which is the default the replicates are not labeled
    plot_min <- ggplot(subset(final_data, TargetRow == min_diff_row_index), aes(x = comparison, y = AverageIntensity, fill = comparison)) +
    geom_violin(trim = FALSE) +
    geom_signif(comparisons = list(c("Ezh2_cKO", "WT")), annotations = paste("p =", formatC(p_values_min, format = "e", digits = 2)), y_position = y_max - 15) +
    geom_jitter(data = subset(bio_rep_average, TargetRow == min_diff_row_index),  aes(x = comparison, y = AverageIntensity, group = comparison), width = 0.2, size = 1, colour = "black") +
    labs(x = "Genotype",
         y = paste("Normalized", sel_frame, "\nmean intensity")) +  # Dynamic y-axis label
    ylim(-10, 150) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = title_min) +
    theme_minimal() + 
    centered_title_theme +
    theme(
      legend.position = "none")
  
    plot_max <- ggplot(subset(final_data, TargetRow == max_diff_row_index), aes(x = comparison, y = AverageIntensity, fill = comparison)) +
    geom_violin(trim = FALSE) +
    geom_signif(comparisons = list(c("Ezh2_cKO", "WT")), annotations = paste("p =", formatC(p_values_max, format = "e", digits = 2)), y_position = y_max - 15) +
    geom_jitter(data = subset(bio_rep_average, TargetRow == max_diff_row_index),  aes(x = comparison, y = AverageIntensity, group = comparison), width = 0.2, size = 1, colour = "black") +
    labs(x = "Genotype", 
         y = paste("Normalized", sel_frame, "\nmean intensity")) +  # Dynamic y-axis label
    scale_fill_brewer(palette = "Set1") +
    ylim(-10, 150) +
    labs(title = title_max) +
    theme_minimal() +
    centered_title_theme +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
} else {
  #TODO add a case where each image is treated individually 
  
  
}





# Combine the plots side by side
combined_plot <- plot_min + plot_max

# Use align_patches to align the plots and ensure equal plot areas
violin_max_min <- combined_plot + plot_layout(ncol = 2)


combined_plot <- frame_line / avg_difference_plot / violin_max_min

# Adjust layout with specified heights for each plot
combined_plot <- combined_plot + plot_layout(heights = c(2.2, 1, 2))

# Display the combined plot
return(combined_plot)

}
