# ============================================================================
# FOREST PLOT GENERATOR FROM CSV FILES
# Standalone script to create publication-ready forest plots
# ============================================================================

# ── LIBRARIES ────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(viridis)
library(purrr)
library(readr)

# ── COLORBLIND-FRIENDLY LAG COLORS ───────────────────────────────────────────
# Wong 2011 palette - highly distinguishable, colorblind-safe
lag_colors <- c(
  "1" = "#E69F00",  # Orange
  "2" = "#56B4E9",  # Sky Blue
  "3" = "#009E73",  # Bluish Green
  "4" = "#F0E442",  # Yellow
  "5" = "#D55E00"   # Vermillion
)

# ── CONFIGURATION ────────────────────────────────────────────────────────────
# Directory containing CSV files
csv_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/Datasets/gentoo-abundance-model/forest_plot_data"

results_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/results"

# Output directory for plots
output_dir <- file.path(results_dir, "Forest_Plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# HELPER FUNCTION: Order home range sizes numerically
# ============================================================================
order_home_range <- function(hr_vector) {
  numeric_values <- as.numeric(gsub("km", "", hr_vector))
  factor(hr_vector, levels = unique(hr_vector[order(numeric_values)]))
}

# ============================================================================
# LOAD AND COMBINE ALL CSV FILES
# ============================================================================

cat("\n", rep("=", 70), "\n")
cat("FOREST PLOT GENERATOR\n")
cat(rep("=", 70), "\n\n")

cat("Reading CSV files from:", csv_dir, "\n")

# Get all CSV files
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No CSV files found in directory: ", csv_dir)
}

cat(sprintf("Found %d CSV files\n", length(csv_files)))

# Read and combine all CSVs
all_data <- csv_files %>%
  set_names(basename(.)) %>%
  map_dfr(~ {
    cat(sprintf("  Reading: %s\n", basename(.x)))
    read_csv(.x, show_col_types = FALSE)
  }, .id = "source_file")

cat(sprintf("\nTotal rows loaded: %d\n", nrow(all_data)))

# ============================================================================
# PREPARE DATA FOR PLOTTING
# ============================================================================

cat("\nPreparing data for plotting...\n")

# Standardize column names (handle variations)
plot_data <- all_data %>%
  rename_with(~ gsub(" ", "", .x)) %>%  # Remove spaces
  rename_with(tolower) %>%               # Lowercase
  mutate(
    # Standardize column names
    metric = if("metric" %in% names(.)) metric else NA_character_,
    homerange = if("homerange" %in% names(.)) homerange 
    else if("homerangesize" %in% names(.)) homerangesize
    else if("homerrange" %in% names(.)) homerrange
    else NA_character_,
    threshold = if("threshold" %in% names(.)) threshold else NA_real_,
    lag = if("lag" %in% names(.)) lag else NA_integer_,
    coefficient = if("coefficient" %in% names(.)) coefficient else NA_real_,
    stderr = if("stderr" %in% names(.)) stderr 
    else if("stderror" %in% names(.)) stderror
    else NA_real_,
    cilower = if("cilower" %in% names(.)) cilower 
    else if("ci_lower" %in% names(.)) ci_lower
    else coefficient - 1.96 * stderr,
    ciupper = if("ciupper" %in% names(.)) ciupper 
    else if("ci_upper" %in% names(.)) ci_upper
    else coefficient + 1.96 * stderr
  ) %>%
  # Clean up home range format
  mutate(
    homerange = as.character(homerange),
    homerange = if_else(str_detect(homerange, "km$"), homerange, paste0(homerange, "km"))
  ) %>%
  # Order home ranges properly
  mutate(
    hr_size = order_home_range(homerange),
    lag = factor(lag, levels = as.character(1:5)),
    threshold = as.character(threshold)
  ) %>%
  # Classify metrics
  mutate(
    metric_type = if_else(
      str_detect(metric, "Variability|SD|sd_|_SD"), 
      "Variability", 
      "Main"
    )
  )

# ============================================================================
# FIX VARIABILITY METRIC NAMES
# ============================================================================

plot_data <- plot_data %>%
  mutate(
    metric = case_when(
      metric == "Variability" & str_detect(source_file, "Concentration") ~ "Sea Ice Concentration Variability",
      metric == "Variability" & str_detect(source_file, "Extent") ~ "Sea Ice Extent Variability",
      metric == "Variability" & str_detect(source_file, "Duration") ~ "Duration Variability",
      metric == "Variability" & str_detect(source_file, "OpenWater") ~ "Open Water Frequency Variability",
      TRUE ~ metric
    )
  )

# Update metric_type classification after renaming
plot_data <- plot_data %>%
  mutate(
    metric_type = if_else(
      str_detect(metric, "Variability"), 
      "Variability", 
      "Main"
    )
  )

cat(sprintf("Prepared %d rows for plotting\n", nrow(plot_data)))
cat(sprintf("Unique metrics: %d\n", n_distinct(plot_data$metric)))
cat(sprintf("Home range sizes: %s\n", paste(levels(plot_data$hr_size), collapse = ", ")))
cat(sprintf("Lags: %s\n", paste(sort(unique(plot_data$lag)), collapse = ", ")))

# ============================================================================
# PATCH: Figures 06 & 07 - Specific Metric Selection and Ordering
# ============================================================================

# Define metric names and order for figures
fig06_metrics <- c(
  "Sea Ice Concentration",
  "Sea Ice Extent", 
  "Duration",
  "Open Water Frequency"
)

fig07_metrics <- c(
  "Sea Ice Concentration Variability",
  "Sea Ice Extent Variability",
  "Duration Variability", 
  "Open Water Frequency Variability"
)

# Filter and prepare data for Figure 06
main_df <- plot_data %>%
  filter(metric_type == "Main") %>%
  filter(metric %in% fig06_metrics) %>%
  filter(threshold %in% c("0.15", "0.3", "0.5")) %>%
  filter(!is.na(coefficient)) %>%
  mutate(
    metric = factor(metric, levels = fig06_metrics),
    threshold = factor(threshold, levels = c("0.15", "0.3", "0.5"))
  )

# Filter and prepare data for Figure 07
var_df <- plot_data %>%
  filter(metric_type == "Variability") %>%
  filter(metric %in% fig07_metrics) %>%
  filter(threshold %in% c("0.15", "0.3", "0.5")) %>%
  filter(!is.na(coefficient)) %>%
  mutate(
    metric = factor(metric, levels = fig07_metrics),
    threshold = factor(threshold, levels = c("0.15", "0.3", "0.5"))
  )



# ── FIGURE 06: Main Metrics (4 rows × 3 columns) ─────────────────────────

if (nrow(main_df) > 0) {
  
  fig06_forest <- ggplot(main_df, aes(x = coefficient, y = hr_size, color = lag)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(
      aes(xmin = cilower, xmax = ciupper),
      height = 0.2, 
      position = position_dodge(width = 0.7),
      na.rm = TRUE
    ) +
    geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
    facet_grid(metric ~ threshold, scales = "free_y") +
    scale_color_manual(name = "Lag (years)", values = lag_colors) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.title.x = element_text(margin = margin(t = 12)),
      axis.title.y = element_text(margin = margin(r = 12))
    ) +
    labs(
      x = "Coefficient (95% CI)",
      y = "Home Range Size (km)"
    )
  
  print(fig06_forest)
  
  ggsave(file.path(output_dir, "Figure_06_Main_Metrics.png"), 
         fig06_forest, width = 12, height = 10, dpi = 300)
  ggsave(file.path(output_dir, "Figure_06_Main_Metrics.pdf"), 
         fig06_forest, width = 12, height = 10)
  
  cat("  Saved: Figure_06_Main_Metrics.png/pdf\n")
}

# ── FIGURE 07: Variability Metrics (4 rows × 3 columns) ──────────────────

if (nrow(var_df) > 0) {
  
  # Strip "Variability" from display names for cleaner labels
  var_df_plot <- var_df %>%
    mutate(metric_display = str_remove(metric, " Variability$"))
  
  fig07_forest <- ggplot(var_df_plot, aes(x = coefficient, y = hr_size, color = lag)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(
      aes(xmin = cilower, xmax = ciupper),
      height = 0.2, 
      position = position_dodge(width = 0.7),
      na.rm = TRUE
    ) +
    geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
    facet_grid(metric_display ~ threshold, scales = "free_y") +
    scale_color_manual(name = "Lag (years)", values = lag_colors) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.title.x = element_text(margin = margin(t = 12)),
      axis.title.y = element_text(margin = margin(r = 12))
    ) +
    labs(
      x = "Coefficient (95% CI)",
      y = "Home Range Size (km)"
    )
  
  print(fig07_forest)
  
  ggsave(file.path(output_dir, "Figure_07_Variability_Metrics.png"), 
         fig07_forest, width = 12, height = 10, dpi = 300)
  ggsave(file.path(output_dir, "Figure_07_Variability_Metrics.pdf"), 
         fig07_forest, width = 12, height = 10)
  
  cat("  Saved: Figure_07_Variability_Metrics.png/pdf\n")
}
# ============================================================================
# SUMMARY
# ============================================================================

cat("\n", rep("=", 70), "\n")
cat("FOREST PLOT GENERATION COMPLETE!\n")
cat(rep("=", 70), "\n\n")

cat("Summary:\n")
cat(sprintf("  Total CSV files processed: %d\n", length(csv_files)))
cat(sprintf("  Total data rows: %d\n", nrow(plot_data)))
cat(sprintf("  Unique metrics: %d\n", n_distinct(plot_data$metric)))
cat(sprintf("  Main metrics plots: %d\n", if(exists("forest_main")) 1 else 0))
cat(sprintf("  Variability metrics plots: %d\n", if(exists("forest_var")) 1 else 0))

cat("\nAll plots saved to:", output_dir, "\n")
cat("\nPlot files:\n")
cat("  - forest_main_metrics.png/pdf\n")
cat("  - forest_variability_metrics.png/pdf\n")
cat("  - forest_combined_all_metrics.png/pdf\n")
cat("  - forest_[metric_name].png/pdf (individual)\n")

cat("\n", rep("=", 70), "\n")