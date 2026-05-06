## Generates figures for the example_data/ images and saves them to figures/
## Run from the repo root: Rscript generate_example_figures.R

library(EBImage)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(zoo)
library(RColorBrewer)

source("align_folia_images_connect.R")
source("classify_connectivity_segments.R")
source("elg_igl_normalize.R")

dir.create("figures", showWarnings = FALSE)

pixel_width <- 0.5119  # µm/pixel (20x objective)

# ── colour scheme ────────────────────────────────────────────────────────────
genotype_colors <- brewer.pal(3, "Set1")[1:2]
names(genotype_colors) <- c("WT", "Ezh2 cKO")

channel_colors <- c(DAPI = "#4477AA", p27 = "#EE6677", NeuN = "#228833")


# ── Figure 1: raw images (p27 channel, all three channels side by side) ──────
message("Fig 1: raw images")

img_wt <- readImage("example_data/2018_05_22_s2_3_p27-0025_fused_crop_3_4.tif")
img_ko <- readImage("example_data/2018_05_22_s3_5_p27-0006_fused_crop_3_4.tif")

# normalise each channel 0-1 for display
norm_frame <- function(img, ch) {
  x <- imageData(img)[,,ch]
  (x - min(x)) / (max(x) - min(x))
}

make_image_df <- function(img, label, channels = c("p27", "NeuN", "DAPI")) {
  lapply(seq_along(channels), function(i) {
    m <- norm_frame(img, i)
    expand.grid(x = seq_len(ncol(m)), y = seq_len(nrow(m))) |>
      mutate(value = as.vector(t(m)),
             channel = channels[i],
             genotype = label)
  }) |> bind_rows()
}

img_df <- bind_rows(
  make_image_df(img_wt, "WT"),
  make_image_df(img_ko, "Ezh2 cKO")
) |>
  mutate(channel  = factor(channel,  levels = c("p27","NeuN","DAPI")),
         genotype = factor(genotype, levels = c("WT","Ezh2 cKO")))

fig1 <- ggplot(img_df, aes(x = x, y = -y, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colours = c("black","white"), guide = "none") +
  facet_grid(genotype ~ channel) +
  coord_fixed() +
  labs(title = "Example images — P7 mouse cerebellum",
       subtitle = "Pial surface at top; deeper layers toward bottom",
       x = NULL, y = NULL) +
  theme_void() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey40"),
    plot.margin   = margin(8, 8, 8, 8)
  )

ggsave("figures/fig1_raw_images.png", fig1, width = 7, height = 5, dpi = 150)
message("  saved figures/fig1_raw_images.png")


# ── load & process both images ───────────────────────────────────────────────
message("Loading images through pipeline ...")

aligned <- align_folia_images_connect(
  "example_data/",
  frame_names          = c("p27", "NeuN", "DAPI"),
  reference_frame      = "DAPI",
  pixel_segment_length = 150
)

# position 15 of the filename encodes the mouse number;
# mouse 3 = WT, mouse 5 = Ezh2 cKO (from the study genotyping)
aligned <- aligned |>
  mutate(
    mouse_num = substr(filename, 15, 15),
    genotype  = case_when(mouse_num == "3" ~ "WT",
                          mouse_num == "5" ~ "Ezh2 cKO",
                          TRUE ~ "unknown"),
    mouse_id  = paste0("M", mouse_num)
  )


# ── Figure 2: raw row-mean profiles ──────────────────────────────────────────
message("Fig 2: raw profiles")

fig2 <- aligned |>
  filter(MaxCorrelation > 0.5) |>
  mutate(channel = factor(Frame, levels = c("p27","NeuN","DAPI"))) |>
  ggplot(aes(x = Row, y = RowMean,
             group = paste0(filename, Frame),
             colour = genotype)) +
  geom_line(alpha = 0.6, linewidth = 0.5) +
  scale_colour_manual(values = genotype_colors) +
  facet_wrap(~channel, ncol = 1, scales = "free_y") +
  labs(title = "Raw row-mean intensity profiles",
       x     = "Row (pixels from image top)",
       y     = "Mean pixel intensity",
       colour = "Genotype") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave("figures/fig2_raw_profiles.png", fig2, width = 5, height = 7, dpi = 150)
message("  saved figures/fig2_raw_profiles.png")


# ── pipeline: filter → prune → shift → split → normalize ────────────────────
message("Running pipeline ...")

filtered <- aligned |> filter(MaxCorrelation > 0.5)

row_coverage <- filtered |>
  group_by(Row) |>
  summarise(rel = n_distinct(filename) / max(n_distinct(filename)), .groups = "drop") |>
  filter(rel >= 0.5)   # relaxed for 2-file example

pruned <- inner_join(filtered, row_coverage, by = "Row")

avg_p27 <- pruned |>
  filter(Frame == "p27") |>
  group_by(Row) |>
  summarise(AvgRowMean = mean(RowMean, na.rm = TRUE), .groups = "drop")

find_local_minimum <- function(df) {
  d   <- diff(df$AvgRowMean)
  cps <- which(d[-1] > 0 & d[-length(d)] < 0)
  cps[which.min(abs(cps - nrow(df) / 2))]
}

origin_row <- find_local_minimum(avg_p27)

shifted <- pruned |>
  mutate(Row_shift = Row - origin_row,
         Row_micron = Row_shift * pixel_width)

neg_half <- shifted |> filter(Row_shift <= 0) |>
  mutate(Row_shift = abs(Row_shift), xpslit = "negative")
pos_half <- shifted |> filter(Row_shift >= 0) |>
  mutate(xpslit = "positive")

split_data <- bind_rows(neg_half, pos_half) |>
  rename(Row_shift_scale = Row_shift) |>
  group_by(Frame, filename, xpslit) |>
  mutate(
    RescaledIntensity = (RowMean - quantile(RowMean, 0.05, na.rm = TRUE)) /
      (quantile(RowMean, 0.95, na.rm = TRUE) -
         quantile(RowMean, 0.05, na.rm = TRUE)) * 100
  ) |>
  ungroup()

normalized <- egl_igl_normalize(split_data, DAPI_max_site = "igl", smooth = TRUE)


# ── Figure 3: normalized profiles with average overlay ───────────────────────
message("Fig 3: normalized profiles")

avg_norm <- normalized |>
  group_by(Row_shift_scale, Frame, genotype) |>
  summarise(avg = mean(RescaledIntensityInRange, na.rm = TRUE), .groups = "drop") |>
  mutate(Row_micron = Row_shift_scale * pixel_width,
         channel    = factor(Frame, levels = c("p27","NeuN","DAPI")))

strip_norm <- normalized |>
  mutate(Row_micron = Row_shift_scale * pixel_width,
         channel    = factor(Frame, levels = c("p27","NeuN","DAPI")))

fig3 <- ggplot() +
  geom_line(data = strip_norm,
            aes(x = Row_micron, y = RescaledIntensityInRange,
                group = paste0(filename, xpslit),
                colour = genotype),
            alpha = 0.2, linewidth = 0.4) +
  geom_line(data = avg_norm,
            aes(x = Row_micron, y = avg,
                colour = genotype, group = genotype),
            linewidth = 1.2) +
  scale_colour_manual(values = genotype_colors) +
  facet_wrap(~channel, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 200)) +
  labs(title    = "EGL/IGL-normalized intensity profiles",
       subtitle = "Bold line = average; faint lines = individual segments",
       x        = "Distance from p27 minimum (µm)",
       y        = "Normalized intensity",
       colour   = "Genotype") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave("figures/fig3_normalized_profiles.png", fig3, width = 5, height = 7, dpi = 150)
message("  saved figures/fig3_normalized_profiles.png")


# ── Figure 4: connectivity metrics ───────────────────────────────────────────
# conn_small_num is the density of small connected components per row —
# a proxy for individual nucleus count. Useful for comparing nuclear density
# across layers and between genotypes.
message("Fig 4: connectivity metrics")

conn_data <- normalized |>
  filter(Frame == "DAPI") |>
  mutate(Row_micron = Row_shift_scale * pixel_width) |>
  group_by(Row_micron, genotype) |>
  summarise(
    avg_conn_small_num = mean(conn_small_num, na.rm = TRUE),
    avg_conn_large     = mean(conn_large,     na.rm = TRUE),
    avg_thr_mean       = mean(thr_mean,       na.rm = TRUE),
    .groups = "drop"
  )

p4a <- ggplot(conn_data, aes(x = Row_micron, y = avg_conn_small_num, colour = genotype)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = genotype_colors) +
  coord_cartesian(xlim = c(0, 200)) +
  labs(y = "Small nucleus density", x = NULL, colour = "Genotype") +
  theme_minimal(base_size = 10) + theme(legend.position = "none")

p4b <- ggplot(conn_data, aes(x = Row_micron, y = avg_thr_mean, colour = genotype)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = genotype_colors) +
  coord_cartesian(xlim = c(0, 200)) +
  labs(y = "DAPI threshold fraction", x = "Distance from p27 minimum (µm)",
       colour = "Genotype") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")

fig4 <- p4a / p4b +
  plot_annotation(
    title    = "DAPI connectivity metrics across cortical depth",
    subtitle = "Derived from connected-component labelling of thresholded images"
  )

ggsave("figures/fig4_connectivity_metrics.png", fig4, width = 6, height = 5, dpi = 150)
message("  saved figures/fig4_connectivity_metrics.png")


# ── Figure 5: p27 overlay — all three normalizations side by side ─────────────
message("Fig 5: p27 normalization comparison")

p27_raw <- normalized |>
  filter(Frame == "p27") |>
  mutate(Row_micron = Row_shift_scale * pixel_width)

avg_raw  <- p27_raw |> group_by(Row_micron, genotype) |>
  summarise(avg = mean(RowMean,               na.rm = TRUE), .groups = "drop")
avg_resc <- p27_raw |> group_by(Row_micron, genotype) |>
  summarise(avg = mean(RescaledIntensity,     na.rm = TRUE), .groups = "drop")
avg_norm <- p27_raw |> group_by(Row_micron, genotype) |>
  summarise(avg = mean(RescaledIntensityInRange, na.rm = TRUE), .groups = "drop")

make_panel <- function(strip_y, avg_y, ylabel, title) {
  ggplot() +
    geom_line(data = p27_raw,
              aes(x = Row_micron, y = .data[[strip_y]],
                  group = paste0(filename, xpslit), colour = genotype),
              alpha = 0.15, linewidth = 0.4) +
    geom_line(data = if (avg_y == "RowMean") avg_raw
              else if (avg_y == "RescaledIntensity") avg_resc
              else avg_norm,
              aes(x = Row_micron, y = avg, colour = genotype, group = genotype),
              linewidth = 1.3) +
    scale_colour_manual(values = genotype_colors) +
    coord_cartesian(xlim = c(0, 200)) +
    labs(title = title, x = "Distance from p27 minimum (µm)",
         y = ylabel, colour = "Genotype") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", plot.title = element_text(size = 9))
}

fig5 <- make_panel("RowMean",               "RowMean",               "Raw intensity",        "1. Raw") +
        make_panel("RescaledIntensity",     "RescaledIntensity",     "Per-strip rescaled",   "2. Per-strip (5–95%)") +
        make_panel("RescaledIntensityInRange", "RescaledIntensityInRange", "EGL/IGL normalized", "3. EGL/IGL anchored") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  plot_annotation(
    title   = "p27 — effect of normalization strategy",
    caption = "Single image per group; for illustration only"
  )

ggsave("figures/fig5_p27_normalization_comparison.png", fig5, width = 10, height = 4, dpi = 150)
message("  saved figures/fig5_p27_normalization_comparison.png")

message("\nAll figures saved to figures/")
