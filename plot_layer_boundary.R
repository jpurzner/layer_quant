generate_combined_plot <- function(data, extra_anno = NULL, sel_frame = "p27", window_size = 5) {

# Load the patchwork package
library(patchwork)
library(lmerTest)
library(RColorBrewer)
library(ggplot2)
library(dplyr)


# set variables for colour   
genotype_colors <- brewer.pal(2, "Set1")
color_mapping <- setNames(genotype_colors, c("WT", "Ezh2 cKO"))
dashed_color <- brewer.pal(4, "Set1")[4]

# First, calculate the average for each Row_shift_scale within each facet
averages <- data %>%
  dplyr::ungroup() %>%
  dplyr::filter(Frame == sel_frame) %>%
  #dplyr::mutate(Row_shift_scale_micron = Row_shift_scale* pixel_width) %>%
  dplyr::group_by(Row_shift_scale_micron, genotype, Frame ) %>%
  dplyr::summarise(average_intensity = mean(RescaledIntensity, na.rm = FALSE)) %>%
  dplyr::mutate(genotype = factor(genotype, levels = c("WT", "Ezh2_cKO")),
                genotype = recode(genotype, "Ezh2_cKO" = "Ezh2 cKO")) 

averages_wide <- averages %>%
  pivot_wider(names_from = genotype, values_from = average_intensity) %>%
  dplyr::filter(Row_shift_scale_micron > 4 & Row_shift_scale_micron < 200) %>%
  dplyr::mutate(Difference = `Ezh2 cKO` - `WT`) # Calculate the difference between genotypes

max_diff_index <- which.max(averages_wide$Difference)

min_diff_index <- which.min(averages_wide$Difference)

# Retrieve the corresponding Row_shift_scale_micron
max_diff_row_scale <- averages_wide$Row_shift_scale_micron[max_diff_index]
min_diff_row_scale <- averages_wide$Row_shift_scale_micron[min_diff_index]



# generate the line plot 
frame_line <-  data %>%
  dplyr::mutate(genotype = factor(genotype, levels = c("WT", "Ezh2_cKO")),
                genotype = recode(genotype, "Ezh2_cKO" = "Ezh2 cKO")) %>%
  dplyr::filter(Frame == sel_frame) %>%
  ggplot() + 
  geom_line(aes(x = Row_shift_scale_micron, y = RescaledIntensity, group = paste0(filename, xpslit)), 
            colour = "black", alpha = 0.05) +
  geom_line(data = averages, aes(x = Row_shift_scale_micron, y = average_intensity, group = genotype, color = genotype), 
            size = 1) +
  geom_vline(xintercept = max_diff_row_scale, linetype = "dashed", color = dashed_color) +
  geom_vline(xintercept = min_diff_row_scale, linetype = "dashed", color = dashed_color) +
  scale_color_manual(values = color_mapping) +
  labs(title = "",
       x = "",
       y = "Normalized p27\nmean intensity") +
  xlim(0, 150) +
  ylim(0,150) +
  facet_grid(genotype ~ Frame) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = -1, unit = "cm")  # Negative bottom margin to bring plots closer
  )


# Plot the difference
avg_difference_plot <- ggplot(averages_wide, aes(x = Row_shift_scale_micron, y = Difference)) +
  geom_line() +
  geom_vline(xintercept = max_diff_row_scale, linetype = "dashed", color = dashed_color) +
  geom_vline(xintercept = min_diff_row_scale, linetype = "dashed", color = dashed_color) +
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


target_rows <- c( max_diff_index, min_diff_index)

# find the micron values nearest to the target rows

expanded_ranges <- lapply(target_rows, function(x) {
  data.frame(Row_shift_scale = (x- window_size):(x + window_size), TargetRow = x)
}) %>% bind_rows()

# Merge expanded_ranges with the filtered data
filtered_expanded_data <- data %>%
  filter(Frame == sel_frame) %>%
  inner_join(expanded_ranges, by = "Row_shift_scale")


averaged_data <- filtered_expanded_data %>%
  dplyr::group_by(filename, xpslit, genotype, TargetRow) %>%
  dplyr::summarize(AverageIntensity = mean(RescaledIntensity, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Merge with the original data to get the associated metadata
merged_data <- dplyr::left_join(averaged_data, 
                                data_pruned_shift_rescaled_split_all_micron %>% 
                                  dplyr::select(filename, xpslit, genotype, Row_shift_scale, dataset, Row_shift_scale_micron), 
                                by = c("filename", "xpslit", "genotype", "TargetRow" = "Row_shift_scale"))

# Ensure unique rows by selecting the first occurrence of each group
final_data <- merged_data %>%
  dplyr::group_by(filename, xpslit, genotype, TargetRow) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

final_data_anno <- dplyr::left_join(final_data, extra_anno, by = c("filename", "xpslit"))

final_data_anno$dataset <- as.factor(final_data_anno$dataset)
final_data_anno$CerebellumPart <- as.factor(final_data_anno$CerebellumPart)
final_data_anno$dataset <- as.factor(final_data_anno$dataset)

mouse_average <- final_data_anno %>%
  dplyr::group_by(Mouse, dataset, genotype, TargetRow, Row_shift_scale_micron) %>%
  dplyr::summarize(AverageIntensity = mean(AverageIntensity, na.rm = TRUE)) %>%
  dplyr::ungroup()

mouse_average$dataset <- as.factor(mouse_average$dataset)
#mouse_average$CerebellumPart <- as.factor(mouse_average$CerebellumPart)
mouse_average$dataset <- as.factor(mouse_average$dataset)

max_diff_final_data_anno <- final_data_anno %>%
  dplyr::filter(TargetRow == max_diff_index)                               

model <- lmer(AverageIntensity ~ genotype +  ( CerebellumPart | Mouse ), data = max_diff_final_data_anno)

summary(model)
model_summary <- summary(model)

p_values_max <- model_summary$coefficients[, "Pr(>|t|)"]


max_diff_final_data_anno <- final_data_anno %>%
  dplyr::filter(TargetRow == min_diff_index)                               

model <- lmer(AverageIntensity ~ genotype +  ( CerebellumPart | Mouse ), data = max_diff_final_data_anno)

summary(model)
model_summary <- summary(model)

p_values_min <- model_summary$coefficients[, "Pr(>|t|)"]

micron_min <- unique(final_data[final_data$TargetRow == min_diff_index, "Row_shift_scale_micron"])
micron_max <- unique(final_data[final_data$TargetRow == max_diff_index, "Row_shift_scale_micron"])

# Create custom titles
title_min <- paste("Micrometers from\npial boundary:", round(micron_min, 2))
title_max <- paste("Micrometers from\npial boundary:", round(micron_max, 2))

centered_title_theme <- theme(plot.title = element_text(hjust = 0.5, size = 12))  # Adjust the size as needed

# Create the two plots as before
plot_min <- ggplot(subset(final_data, TargetRow == min_diff_index), aes(x = genotype, y = AverageIntensity, fill = genotype)) +
  geom_violin(trim = FALSE) +
  geom_signif(comparisons = list(c("Ezh2_cKO", "WT")), annotations = paste("p =", formatC(p_values_min["genotypeWT"], format = "e", digits = 2)), y_position = max(final_data$AverageIntensity) * 0.9) +
  labs(x = "Genotype", y = "Normalized p27\nmean intensity") +
  ylim(-10, 150) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = title_35) +
  theme_minimal() + 
  centered_title_theme +
  theme(
    legend.position = "none")

plot_max <- ggplot(subset(final_data, TargetRow == max_diff_index), aes(x = genotype, y = AverageIntensity, fill = genotype)) +
  geom_violin(trim = FALSE) +
  geom_signif(comparisons = list(c("Ezh2_cKO", "WT")), annotations = paste("p =", formatC(p_values_max["genotypeWT"], format = "e", digits = 2)), y_position = max(final_data$AverageIntensity) * 0.9) +
  labs(x = "Genotype", y = "Normalized p27\nmean intensity") +
  scale_fill_brewer(palette = "Set1") +
  ylim(-10, 150) +
  labs(title = title_82) +
  theme_minimal() +
  centered_title_theme +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

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
