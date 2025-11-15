# ============================================================================
# FOREST PLOT GENERATOR FROM CSV FILES from Revision 01 Results
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

cat(sprintf("Prepared %d rows for plotting\n", nrow(plot_data)))
cat(sprintf("Unique metrics: %d\n", n_distinct(plot_data$metric)))
cat(sprintf("Home range sizes: %s\n", paste(levels(plot_data$hr_size), collapse = ", ")))
cat(sprintf("Lags: %s\n", paste(sort(unique(plot_data$lag)), collapse = ", ")))

# ============================================================================
# FOREST PLOT 1: MAIN METRICS
# ============================================================================

cat("\n", rep("=", 70), "\n")
cat("CREATING MAIN METRICS FOREST PLOT\n")
cat(rep("=", 70), "\n\n")

main_df <- plot_data %>%
  filter(metric_type == "Main") %>%
  filter(!is.na(coefficient))

if (nrow(main_df) > 0) {
  
  # Check if thresholds exist
  has_threshold <- length(unique(main_df$threshold)) > 1 || 
    !all(is.na(main_df$threshold))
  
  forest_main <- ggplot(main_df, aes(x = coefficient, y = hr_size, color = lag)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(
      aes(xmin = cilower, xmax = ciupper),
      height = 0.2, 
      position = position_dodge(width = 0.7),
      na.rm = TRUE
    ) +
    geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
    scale_color_viridis_d(name = "Lag (years)", end = 0.8) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) +
    labs(
      title = "Forest Plot: Main Sea Ice Metrics",
      x = "Coefficient (95% CI)",
      y = "Home Range Size (km)"
    )
  
  # Add faceting
  if (has_threshold) {
    forest_main <- forest_main +
      facet_grid(metric ~ threshold, scales = "free_y")
  } else {
    forest_main <- forest_main +
      facet_wrap(~ metric, scales = "free_y", ncol = 1)
  }
  
  print(forest_main)
  
  ggsave(file.path(output_dir, "forest_main_metrics.png"), 
         forest_main, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "forest_main_metrics.pdf"), 
         forest_main, width = 10, height = 8)
  
  cat("  Saved: forest_main_metrics.png/pdf\n")
  
} else {
  cat("  No main metrics data to plot\n")
}

# ============================================================================
# FOREST PLOT 2: VARIABILITY METRICS
# ============================================================================

cat("\n", rep("=", 70), "\n")
cat("CREATING VARIABILITY METRICS FOREST PLOT\n")
cat(rep("=", 70), "\n\n")

var_df <- plot_data %>%
  filter(metric_type == "Variability") %>%
  filter(!is.na(coefficient))

if (nrow(var_df) > 0) {
  
  has_threshold <- length(unique(var_df$threshold)) > 1 || 
    !all(is.na(var_df$threshold))
  
  forest_var <- ggplot(var_df, aes(x = coefficient, y = hr_size, color = lag)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(
      aes(xmin = cilower, xmax = ciupper),
      height = 0.2, 
      position = position_dodge(width = 0.7),
      na.rm = TRUE
    ) +
    geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
    scale_color_viridis_d(name = "Lag (years)", end = 0.8) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) +
    labs(
      title = "Forest Plot: Variability Metrics",
      x = "Coefficient (95% CI)",
      y = "Home Range Size (km)"
    )
  
  if (has_threshold) {
    forest_var <- forest_var +
      facet_grid(metric ~ threshold, scales = "free_y")
  } else {
    forest_var <- forest_var +
      facet_wrap(~ metric, scales = "free_y", ncol = 1)
  }
  
  print(forest_var)
  
  ggsave(file.path(output_dir, "forest_variability_metrics.png"), 
         forest_var, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "forest_variability_metrics.pdf"), 
         forest_var, width = 10, height = 8)
  
  cat("  Saved: forest_variability_metrics.png/pdf\n")
  
} else {
  cat("  No variability metrics data to plot\n")
}

# ============================================================================
# FOREST PLOT 3: COMBINED (ALL METRICS)
# ============================================================================

cat("\n", rep("=", 70), "\n")
cat("CREATING COMBINED FOREST PLOT\n")
cat(rep("=", 70), "\n\n")

combined_df <- plot_data %>%
  filter(!is.na(coefficient))

if (nrow(combined_df) > 0) {
  
  has_threshold <- length(unique(combined_df$threshold)) > 1 || 
    !all(is.na(combined_df$threshold))
  
  forest_combined <- ggplot(combined_df, aes(x = coefficient, y = hr_size, color = lag)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(
      aes(xmin = cilower, xmax = ciupper),
      height = 0.2, 
      position = position_dodge(width = 0.7),
      na.rm = TRUE
    ) +
    geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
    scale_color_viridis_d(name = "Lag (years)", end = 0.8) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    ) +
    labs(
      title = "Forest Plot: All Sea Ice Metrics",
      x = "Coefficient (95% CI)",
      y = "Home Range Size (km)"
    )
  
  if (has_threshold) {
    forest_combined <- forest_combined +
      facet_grid(metric ~ threshold, scales = "free_y")
  } else {
    forest_combined <- forest_combined +
      facet_wrap(~ metric, scales = "free_y", ncol = 2)
  }
  
  print(forest_combined)
  
  ggsave(file.path(output_dir, "forest_combined_all_metrics.png"), 
         forest_combined, width = 14, height = 12, dpi = 300)
  ggsave(file.path(output_dir, "forest_combined_all_metrics.pdf"), 
         forest_combined, width = 14, height = 12)
  
  cat("  Saved: forest_combined_all_metrics.png/pdf\n")
}

# ============================================================================
# FOREST PLOT 4: INDIVIDUAL METRIC PLOTS
# ============================================================================

cat("\n", rep("=", 70), "\n")
cat("CREATING INDIVIDUAL METRIC PLOTS\n")
cat(rep("=", 70), "\n\n")

# Derive base metric names
base_metrics <- plot_data %>%
  pull(metric) %>%
  str_remove(" Variability$") %>%
  str_remove(" SD$") %>%
  str_remove("_SD$") %>%
  str_remove(" sd$") %>%
  str_remove("_sd$") %>%
  unique() %>%
  na.omit()

cat(sprintf("Creating plots for %d base metrics\n", length(base_metrics)))

for (bm in base_metrics) {
  
  # Find all variations of this metric
  df_sub <- plot_data %>%
    filter(
      metric == bm | 
        str_detect(metric, paste0("^", bm, " (Variability|SD|sd)")) |
        str_detect(metric, paste0("^", bm, "_(SD|sd)"))
    ) %>%
    filter(!is.na(coefficient))
  
  if (nrow(df_sub) == 0) next
  
  cat(sprintf("  Creating plot for: %s\n", bm))
  
  has_threshold <- length(unique(df_sub$threshold)) > 1 || 
    !all(is.na(df_sub$threshold))
  
  p <- ggplot(df_sub, aes(x = coefficient, y = hr_size, color = lag)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(
      aes(xmin = cilower, xmax = ciupper),
      height = 0.2, 
      position = position_dodge(width = 0.7),
      na.rm = TRUE
    ) +
    geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
    scale_color_viridis_d(name = "Lag (years)", end = 0.8) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      axis.title.y.left = element_blank(),
      axis.text.y.left = element_text(margin = margin(r = 10)),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0("Forest Plot: ", bm),
      x = "Coefficient (95% CI)",
      y = "Home Range Size (km)"
    )
  
  if (has_threshold) {
    p <- p + facet_grid(metric ~ threshold, scales = "free_y", switch = "y")
  } else {
    p <- p + facet_wrap(~ metric, scales = "free_y", ncol = 1, switch = "y")
  }
  
  print(p)
  
  safe_name <- gsub("[^[:alnum:]_]", "_", bm)
  ggsave(
    filename = file.path(output_dir, paste0("forest_", safe_name, ".png")),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  ggsave(
    filename = file.path(output_dir, paste0("forest_", safe_name, ".pdf")),
    plot = p,
    width = 10,
    height = 6
  )
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
cat(sprintf("  Individual metric plots: %d\n", length(base_metrics)))

cat("\nAll plots saved to:", output_dir, "\n")
cat("\nPlot files:\n")
cat("  - forest_main_metrics.png/pdf\n")
cat("  - forest_variability_metrics.png/pdf\n")
cat("  - forest_combined_all_metrics.png/pdf\n")
cat("  - forest_[metric_name].png/pdf (individual)\n")

cat("\n", rep("=", 70), "\n")
