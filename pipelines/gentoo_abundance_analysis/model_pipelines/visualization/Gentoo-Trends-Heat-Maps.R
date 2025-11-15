# ── Complete Updated Script: Fix .csv in Variability Names & Full Heatmap Grid ──

# 1. Load required libraries
library(tidyverse)   # readr, dplyr, ggplot2, purrr, etc.
library(stringr)     # str_remove, str_detect
library(tidyr)       # complete()
library(viridis)     # scale_fill_viridis()

# 2. Define data directory and read all CSVs
data_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/data/visualization"
csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

df_raw <- csv_files %>%
  set_names() %>%
  map_dfr(read_csv, .id = "source_file")

# 3. Disambiguate “Variability” metrics using filenames, stripping .csv
df2 <- df_raw %>%
  mutate(
    file_base = basename(source_file),
    # remove the .csv extension first
    file_base = str_remove(file_base, "\\.csv$"),
    # derive base_metric from filename for Variability rows
    base_metric = if_else(
      Metric == "Variability",
      # keep text before "_Variability"
      str_remove(file_base, "_Variability.*") %>%
        str_replace_all("_", " "),  # underscores → spaces
      NA_character_
    ),
    # build the final metric name
    metric = if_else(
      Metric == "Variability",
      paste0(base_metric, " Variability"),
      Metric
    )
  ) %>%
  # 4. Clean up and rename columns for plotting
  select(-Metric, -file_base) %>%
  rename(
    hr_size   = `Home Range Size`,
    lag       = Lag,
    threshold = Threshold,
    estimate  = Coefficient
  ) %>%
  # 5. Convert to factors with known levels
  mutate(
    hr_size   = factor(hr_size, levels = sort(unique(hr_size))),
    lag       = factor(lag,     levels = as.character(1:5)),
    threshold = factor(threshold, levels = c(0.15, 0.30, 0.50))
  ) %>%
  # 6. Complete the full grid of metric × threshold × hr_size × lag
  complete(metric, threshold, hr_size, lag)

# 7. Build and display the heatmap
heatmap_plot <- ggplot(df2, aes(x = lag, y = hr_size, fill = estimate)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(
    aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
    size = 3.5, color = "white"
  ) +
  scale_fill_viridis(
    option    = "D",
    direction = -1,
    begin     = 0.15,
    end       = 0.75,
    na.value  = "grey95",
    name      = "Coef."
  ) +
  facet_grid(metric ~ threshold, scales = "free_y") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    strip.text       = element_text(face = "bold", size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y      = element_text(size = 11),
    axis.title       = element_text(size = 13),
    legend.position  = "right",
    legend.background= element_blank(),
    plot.margin      = margin(t = 10, r = 15, b = 10, l = 15)
  ) +
  labs(
    title = "Effect Sizes of Winter Sea Ice Metrics on Gentoo Penguin Growth",
    x     = "Lag (years)",
    y     = "Home‑Range Size (km)"
  )

print(heatmap_plot)

# 8. (Optional) Save to file
ggsave(
  filename = "heatmap_full_grid.png",
  plot     = heatmap_plot,
  width    = 12,
  height   = 8,
  dpi      = 300
)

# ── 1. Main Metrics Heatmap ───────────────────────────────────────────────────

main_df <- df2 %>%
  filter(!str_detect(metric, "Variability"))

main_plot <- ggplot(main_df, aes(x = lag, y = hr_size, fill = estimate)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(
    aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
    size = 3.5, color = "white"
  ) +
  scale_fill_viridis(option = "D", direction = -1,
                     begin = 0.15, end = 0.75,
                     na.value = "grey95", name = "Coef.") +
  facet_grid(metric ~ threshold, scales = "free_y") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    strip.text       = element_text(face = "bold", size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y      = element_text(size = 11),
    legend.position  = "right",
    legend.background= element_blank(),
    plot.margin      = margin(t = 10, r = 15, b = 10, l = 15)
  ) +
  labs(
    title = "Main Sea Ice Metrics Coefficients",
    x     = "Lag (years)",
    y     = "Home‑Range Size (km)"
  )

print(main_plot)
ggsave("heatmap_main_metrics.png", main_plot, width=12, height=8, dpi=300)


# ── 2. Variability‐Only Heatmap ────────────────────────────────────────────────

var_df <- df2 %>%
  filter(str_detect(metric, "Variability"))

var_plot <- ggplot(var_df, aes(x = lag, y = hr_size, fill = estimate)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(
    aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
    size = 3.5, color = "white"
  ) +
  scale_fill_viridis(option = "D", direction = -1,
                     begin = 0.15, end = 0.75,
                     na.value = "grey95", name = "Coef.") +
  facet_grid(metric ~ threshold, scales = "free_y") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    strip.text       = element_text(face = "bold", size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y      = element_text(size = 11),
    legend.position  = "right",
    legend.background= element_blank(),
    plot.margin      = margin(t = 10, r = 15, b = 10, l = 15)
  ) +
  labs(
    title = "Variability Metrics Coefficients",
    x     = "Lag (years)",
    y     = "Home‑Range Size (km)"
  )

print(var_plot)
ggsave("heatmap_variability_metrics.png", var_plot, width=12, height=8, dpi=300)




# ── Revised save_nature_figure: use mm units directly ───────────

save_nature_figure <- function(plot, name, width_mm, height_mm){
  # Vector PDF
  ggsave(
    filename  = paste0(name, ".pdf"),
    plot      = plot,
    device    = cairo_pdf,
    width     = width_mm,
    height    = height_mm,
    units     = "mm",       # mm instead of inches
    dpi       = 300,        # 300 dpi for vector figure
    limitsize = FALSE
  )
  # High‑res TIFF
  ggsave(
    filename    = paste0(name, ".tiff"),
    plot        = plot,
    device      = "tiff",
    width       = width_mm,
    height      = height_mm,
    units       = "mm",
    dpi         = 600,      # 600 dpi for line art
    compression = "lzw",
    limitsize   = FALSE
  )
}

# ── Usage ──────────────────────────────────────────────────────────

# 1) Full heatmap (single‑column width = 174 mm, height = 234 mm)
save_nature_figure(heat_all, "heatmap_full_grid_nature", 174, 234)

# 2) Main metrics (double‑column width = 84 mm, same height)
save_nature_figure(main_plot, "heatmap_main_metrics_nature", 84, 234)

# 3) Variability metrics (double‑column width = 84 mm)
save_nature_figure(var_plot, "heatmap_variability_metrics_nature", 84, 234)


library(magick)

# 1. Read the TIFF file
img <- image_read("heatmap_full_grid_nature.tiff")

# 2. Inspect its dimensions & resolution
info <- image_info(img)
print(info)
