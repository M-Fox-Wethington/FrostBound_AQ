# ============================================================================
# visualization_helpers.R
# Visualization and output generation for penguin-sea ice analysis
#
# Dependencies: Requires data_processing_helpers.R to be sourced
# Required packages: ggplot2, patchwork, dplyr, data.table
# ============================================================================

#' Create publication-quality scatter plot with GLS fit
#'
#' @param penguin_data Data frame with penguin observations
#' @param gls_model Fitted GLS model object
#' @param intercept_ci Confidence interval for intercept [lower, est, upper]
#' @param slope_ci Confidence interval for slope [lower, est, upper]
#' @param metric String name of metric column in penguin_data
#' @param x_label String label for x-axis
#' @param title String plot title
#' @param region_label String region identifier
#' @return ggplot object
create_plot <- function(penguin_data, gls_model, intercept_ci, slope_ci, 
                        metric, x_label, title, region_label) {
  new_data <- data.frame(overwinter_metric = seq(min(penguin_data[[metric]], na.rm = TRUE),
                                                 max(penguin_data[[metric]], na.rm = TRUE),
                                                 length.out = 100))
  
  colnames(new_data) <- metric
  
  predicted_values <- predict(gls_model, new_data)
  
  new_data$fit <- predicted_values
  new_data$upper <- intercept_ci[3] + slope_ci[3] * new_data[[metric]]
  new_data$lower <- intercept_ci[1] + slope_ci[1] * new_data[[metric]]
  
  plot <- ggplot() +
    geom_point(data = penguin_data, aes_string(x = metric, y = "growth_rate"), 
               color = "black", size = 1, alpha = 0.6) +
    geom_line(data = new_data, aes_string(x = metric, y = "fit"), 
              color = "blue", size = 1) +
    geom_ribbon(data = new_data, aes_string(x = metric, ymin = "lower", ymax = "upper"), 
                alpha = 0.2, fill = "blue") +
    labs(title = paste(title, "-", region_label),
         x = x_label,
         y = "Growth Rate") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0.5, size = 10),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    ylim(NA, 1.5)
  
  return(plot)
}

#' Create side-by-side comparison plots for significant regional results
#'
#' @param results_dir String path to results directory
#' @param penguin_data_path String path to penguin data CSV
#' @param metric_files Character vector of metric file paths
#' @param central_cutoff_lat Numeric latitude threshold for regional classification
#' @return NULL (saves plots to disk)
create_regional_comparison_plots <- function(results_dir, 
                                             penguin_data_path,
                                             metric_files,
                                             central_cutoff_lat = -63.2) {
  
  cat("\n", rep("=", 70), "\n")
  cat("CREATING REGIONAL COMPARISON PLOTS\n")
  cat(rep("=", 70), "\n\n")
  
  # Load penguin data
  penguin_full <- fread(penguin_data_path)
  lat_col <- pick_lat_col(penguin_full)
  
  # Add regions
  penguin_full <- penguin_full %>%
    mutate(
      year = 1970 + season - 1,
      .lat_raw = .data[[lat_col]],
      .lat_sgn = if_else(.lat_raw < 0, .lat_raw, -abs(.lat_raw)),
      region = if_else(.lat_sgn <= central_cutoff_lat, "Central_WAP", "Bransfield")
    ) %>%
    filter(growth_rate <= 3, !is.na(growth_rate))
  
  # Check if significant results exist for both regions
  brans_sig_path <- file.path(results_dir, "Bransfield", "model_results_significant.csv")
  central_sig_path <- file.path(results_dir, "Central_WAP", "model_results_significant.csv")
  
  if (!file.exists(brans_sig_path) & !file.exists(central_sig_path)) {
    cat("No significant results found for either region. Skipping plots.\n")
    return(NULL)
  }
  
  # Load significant results from both regions
  sig_results_list <- list()
  
  if (file.exists(brans_sig_path)) {
    sig_results_list[["Bransfield"]] <- read.csv(brans_sig_path)
  }
  
  if (file.exists(central_sig_path)) {
    sig_results_list[["Central_WAP"]] <- read.csv(central_sig_path)
  }
  
  # Create plot directory
  plot_dir <- file.path(results_dir, "Comparison_Plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Get all unique metric/threshold/homerange/lag combinations across both regions
  all_combos <- bind_rows(sig_results_list) %>%
    select(Metric, HomeRangeSize, Threshold, Lag) %>%
    distinct() %>%
    arrange(Metric, Lag)
  
  cat(sprintf("Found %d unique significant metric combinations\n", nrow(all_combos)))
  
  # Load all metric data
  metric_data_list <- lapply(metric_files, function(f) {
    dat <- fread(f)
    dat$source_file <- basename(f)
    return(dat)
  })
  
  # For each combination, create side-by-side plots
  for (i in 1:nrow(all_combos)) {
    combo <- all_combos[i, ]
    
    cat(sprintf("\nPlot %d/%d: %s, HR=%s, Threshold=%s, Lag=%d\n",
                i, nrow(all_combos), combo$Metric, combo$HomeRangeSize, 
                combo$Threshold, combo$Lag))
    
    plot_list <- list()
    
    # Create plot for each region that has this significant result
    for (region_name in names(sig_results_list)) {
      
      # Check if this region has this specific result
      region_result <- sig_results_list[[region_name]] %>%
        filter(Metric == combo$Metric,
               HomeRangeSize == combo$HomeRangeSize,
               Threshold == combo$Threshold,
               Lag == combo$Lag)
      
      if (nrow(region_result) == 0) {
        cat(sprintf("  Skipping %s (not significant)\n", region_name))
        next
      }
      
      # Get the metric column name from the result
      metric_col <- region_result$MetricColumn[1]
      
      # Find the appropriate metric file
      metric_file <- metric_files[grepl(region_result$File[1], metric_files)]
      if (length(metric_file) == 0) next
      
      # Load and process metric data
      metric_dat <- fread(metric_file) %>%
        filter(HomeRangeSize == combo$HomeRangeSize,
               (is.na(Threshold) | Threshold == combo$Threshold))
      
      # Aggregate to yearly
      metric_yearly <- yearly_aggregate_metric(metric_dat, metric_col) %>%
        mutate(year = year + combo$Lag)
      
      # Merge with penguin data for this region
      penguin_region <- penguin_full %>%
        filter(region == region_name) %>%
        left_join(metric_yearly, by = "year") %>%
        filter(!is.na(metric_value), !is.na(growth_rate))
      
      if (nrow(penguin_region) == 0) {
        cat(sprintf("  No data available for %s after merging\n", region_name))
        next
      }
      
      # Check for sufficient variation in metric
      if (sd(penguin_region$metric_value, na.rm = TRUE) == 0) {
        cat(sprintf("  No variation in metric for %s, skipping\n", region_name))
        next
      }
      
      # Refit the model for plotting
      model <- gls(growth_rate ~ metric_value, 
                   correlation = corAR1(form = ~1 | site_id), 
                   data = penguin_region)
      
      # Try to get confidence intervals with error handling
      conf_int <- tryCatch(
        intervals(model),
        error = function(e) {
          warning(sprintf("Could not compute confidence intervals for %s: %s", 
                          region_name, e$message))
          return(NULL)
        }
      )
      
      if (is.null(conf_int)) {
        cat(sprintf("  Skipping plot for %s (confidence interval error)\n", region_name))
        next
      }
      
      intercept_ci <- conf_int$coef[1, ]
      slope_ci <- conf_int$coef[2, ]
      
      # Create prediction data
      new_data <- data.frame(
        metric_value = seq(min(penguin_region$metric_value, na.rm = TRUE),
                           max(penguin_region$metric_value, na.rm = TRUE),
                           length.out = 100)
      )
      
      new_data$fit <- predict(model, new_data)
      new_data$upper <- intercept_ci[3] + slope_ci[3] * new_data$metric_value
      new_data$lower <- intercept_ci[1] + slope_ci[1] * new_data$metric_value
      
      # Create plot
      p <- ggplot() +
        geom_point(data = penguin_region, 
                   aes(x = metric_value, y = growth_rate), 
                   color = "black", size = 1.5, alpha = 0.6) +
        geom_line(data = new_data, 
                  aes(x = metric_value, y = fit), 
                  color = "blue", linewidth = 1) +
        geom_ribbon(data = new_data, 
                    aes(x = metric_value, ymin = lower, ymax = upper), 
                    alpha = 0.2, fill = "blue") +
        labs(title = region_name,
             subtitle = sprintf("N=%d, p=%.4f", 
                                nrow(penguin_region),
                                region_result$Bonferroni_Adjusted_pValue[1]),
             x = combo$Metric,
             y = "Growth Rate") +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          plot.margin = margin(10, 10, 10, 10)
        ) +
        ylim(NA, 1.5)
      
      plot_list[[region_name]] <- p
    }
    
    # If we have plots for both regions, combine them
    if (length(plot_list) > 0) {
      combined_plot <- wrap_plots(plot_list, ncol = length(plot_list)) +
        plot_annotation(
          title = sprintf("%s vs Growth Rate (Lag %d years)", combo$Metric, combo$Lag),
          subtitle = sprintf("HomeRange: %s, Threshold: %s", 
                             combo$HomeRangeSize, combo$Threshold),
          theme = theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5)
          )
        )
      
      # Save plot
      plot_filename <- sprintf("%s_HR%s_Thr%s_Lag%d.png",
                               gsub(" ", "_", combo$Metric),
                               combo$HomeRangeSize,
                               combo$Threshold,
                               combo$Lag)
      
      tryCatch({
        suppressWarnings({
          ggsave(file.path(plot_dir, plot_filename),
                 plot = combined_plot,
                 width = 4 * length(plot_list),
                 height = 5,
                 dpi = 300)
        })
        cat(sprintf("  Saved: %s\n", plot_filename))
      }, error = function(e) {
        warning(sprintf("Could not save plot %s: %s", plot_filename, e$message))
      })
    }
  }
  
  cat("\n", rep("=", 70), "\n")
  cat("PLOT GENERATION COMPLETE\n")
  cat(rep("=", 70), "\n")
}

#' Create regional comparison summary table
#'
#' @param results_dir String path to results directory
#' @return Data frame with comparison statistics
create_regional_comparison_table <- function(results_dir) {
  
  cat("\n", rep("=", 70), "\n")
  cat("CREATING REGIONAL COMPARISON SUMMARY TABLE\n")
  cat(rep("=", 70), "\n\n")
  
  # Load significant results from both regions
  brans_sig_path <- file.path(results_dir, "Bransfield", "model_results_significant.csv")
  central_sig_path <- file.path(results_dir, "Central_WAP", "model_results_significant.csv")
  
  results_list <- list()
  
  if (file.exists(brans_sig_path)) {
    results_list[["Bransfield"]] <- read.csv(brans_sig_path)
  }
  
  if (file.exists(central_sig_path)) {
    results_list[["Central_WAP"]] <- read.csv(central_sig_path)
  }
  
  if (length(results_list) == 0) {
    cat("No significant results found. Cannot create comparison table.\n")
    return(NULL)
  }
  
  # Combine all results
  all_results <- bind_rows(results_list)
  
  # Create comparison summary
  comparison <- all_results %>%
    group_by(Region, Metric, Lag) %>%
    summarise(
      N_Models = n(),
      Mean_Coefficient = mean(Coefficient, na.rm = TRUE),
      SD_Coefficient = sd(Coefficient, na.rm = TRUE),
      Mean_pValue = mean(pValue, na.rm = TRUE),
      Min_pValue = min(pValue, na.rm = TRUE),
      Mean_N_Obs = mean(N_Obs, na.rm = TRUE),
      Mean_N_Colonies = mean(N_Colonies, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Metric, Lag, Region)
  
  # Save comparison table
  write.csv(comparison, 
            file.path(results_dir, "Regional_Comparison_Summary.csv"),
            row.names = FALSE)
  
  cat("Saved summary table: Regional_Comparison_Summary.csv\n")
  
  # Print key comparisons
  cat("\n=== KEY COMPARISONS ===\n")
  
  # For each metric, show if both regions have significant effects
  metrics_in_both <- comparison %>%
    group_by(Metric, Lag) %>%
    filter(n() == 2) %>%  # Both regions present
    ungroup()
  
  if (nrow(metrics_in_both) > 0) {
    cat("\nMetrics significant in BOTH regions:\n")
    for (met in unique(metrics_in_both$Metric)) {
      cat(sprintf("\n%s:\n", met))
      subset_met <- metrics_in_both %>% filter(Metric == met)
      print(subset_met %>% select(Region, Lag, Mean_Coefficient, Min_pValue, Mean_N_Colonies))
    }
  }
  
  # Show metrics only in one region
  metrics_single <- comparison %>%
    group_by(Metric, Lag) %>%
    filter(n() == 1) %>%
    ungroup()
  
  if (nrow(metrics_single) > 0) {
    cat("\n\nMetrics significant in ONLY ONE region:\n")
    print(metrics_single %>% select(Region, Metric, Lag, Mean_Coefficient, Min_pValue))
  }
  
  cat("\n", rep("=", 70), "\n")
  
  return(comparison)
}

#' Create forest plots for regional comparison of effect sizes
#'
#' @param results_dir String path to results directory
#' @return Data frame with all results used for plotting
create_regional_forest_plots <- function(results_dir) {
  
  cat("\n", rep("=", 70), "\n")
  cat("CREATING REGIONAL COMPARISON FOREST PLOTS\n")
  cat(rep("=", 70), "\n\n")
  
  # Load significant results from both regions
  brans_sig_path <- file.path(results_dir, "Bransfield", "model_results_significant.csv")
  central_sig_path <- file.path(results_dir, "Central_WAP", "model_results_significant.csv")
  
  results_list <- list()
  
  if (file.exists(brans_sig_path)) {
    results_list[["Bransfield"]] <- read.csv(brans_sig_path)
  }
  
  if (file.exists(central_sig_path)) {
    results_list[["Central_WAP"]] <- read.csv(central_sig_path)
  }
  
  if (length(results_list) == 0) {
    cat("No significant results found. Cannot create forest plots.\n")
    return(NULL)
  }
  
  # Combine all results
  all_results <- bind_rows(results_list)
  
  # Create forest plot directory
  forest_dir <- file.path(results_dir, "Forest_Plots")
  dir.create(forest_dir, showWarnings = FALSE)
  
  # Recode Lag column for better labels
  all_results$Lag_Label <- recode(all_results$Lag,
                                  `1` = '1 Year Lag',
                                  `2` = '2 Year Lag',
                                  `3` = '3 Year Lag',
                                  `4` = '4 Year Lag',
                                  `5` = '5 Year Lag')
  
  # Recode Metric names for better labels
  all_results$Metric_Label <- recode(all_results$Metric,
                                     'SIC' = 'Sea Ice Concentration',
                                     'SIC_SD' = 'SIC Variability',
                                     'Extent' = 'Sea Ice Extent',
                                     'Duration' = 'Sea Ice Duration',
                                     'Duration_SD' = 'Duration Variability',
                                     'Persistence' = 'Open Water Frequency',
                                     'Persistence_SD' = 'Open Water Freq. Variability')
  
  # ========================================================================
  # PLOT 1: Combined forest plot for all metrics (faceted by metric)
  # ========================================================================
  cat("\nCreating combined forest plot across all metrics...\n")
  
  combined_forest <- ggplot(all_results, 
                            aes(x = HomeRangeSize, y = Coefficient, 
                                color = Region, shape = Region)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = Coefficient - StdError, 
                      ymax = Coefficient + StdError),
                  position = position_dodge(width = 0.5),
                  width = 0.3,
                  linewidth = 0.5) +
    facet_wrap(~ Metric_Label + Lag_Label, scales = 'free', ncol = 5) +
    scale_color_manual(values = c("Bransfield" = "#2166AC", 
                                  "Central_WAP" = "#B2182B")) +
    scale_shape_manual(values = c("Bransfield" = 16, 
                                  "Central_WAP" = 17)) +
    labs(title = "Regional Comparison: Effect of Sea Ice Metrics on Gentoo Growth Rates",
         subtitle = "Bransfield (South Shetland Islands) vs Central WAP (Gerlache Strait/Anvers)",
         x = "Home Range Size",
         y = "Effect Size (Coefficient ± SE)",
         color = "Region",
         shape = "Region") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  # Save combined forest plot
  ggsave(file.path(forest_dir, "Regional_Comparison_Combined_Forest_Plot.png"),
         combined_forest, width = 18, height = 12, dpi = 300)
  
  ggsave(file.path(forest_dir, "Regional_Comparison_Combined_Forest_Plot.pdf"),
         combined_forest, width = 18, height = 12)
  
  cat("  Saved: Regional_Comparison_Combined_Forest_Plot.png/pdf\n")
  
  # ========================================================================
  # PLOT 2: Individual forest plots for each metric
  # ========================================================================
  cat("\nCreating individual forest plots for each metric...\n")
  
  unique_metrics <- unique(all_results$Metric_Label)
  
  for (metric in unique_metrics) {
    
    metric_data <- all_results %>% filter(Metric_Label == metric)
    
    if (nrow(metric_data) == 0) next
    
    cat(sprintf("  Creating plot for: %s\n", metric))
    
    metric_forest <- ggplot(metric_data,
                            aes(x = HomeRangeSize, y = Coefficient,
                                color = Region, shape = Region)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = Coefficient - StdError,
                        ymax = Coefficient + StdError),
                    position = position_dodge(width = 0.5),
                    width = 0.3,
                    linewidth = 0.7) +
      facet_wrap(~ Lag_Label, scales = 'free_x', nrow = 1) +
      scale_color_manual(values = c("Bransfield" = "#2166AC",
                                    "Central_WAP" = "#B2182B")) +
      scale_shape_manual(values = c("Bransfield" = 16,
                                    "Central_WAP" = 17)) +
      labs(title = paste("Regional Comparison:", metric),
           subtitle = "Effect on Gentoo Penguin Growth Rates",
           x = "Home Range Size",
           y = "Effect Size (Coefficient ± SE)",
           color = "Region",
           shape = "Region") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
    
    # Save individual metric plot
    metric_filename <- gsub(" ", "_", tolower(metric))
    ggsave(file.path(forest_dir, paste0("Regional_Forest_", metric_filename, ".png")),
           metric_forest, width = 12, height = 6, dpi = 300)
    
    ggsave(file.path(forest_dir, paste0("Regional_Forest_", metric_filename, ".pdf")),
           metric_forest, width = 12, height = 6)
  }
  
  # ========================================================================
  # PLOT 3: Forest plot showing only metrics significant in BOTH regions
  # ========================================================================
  cat("\nCreating forest plot for metrics significant in both regions...\n")
  
  # Find metrics present in both regions
  metrics_both <- all_results %>%
    group_by(Metric_Label, Lag_Label, HomeRangeSize, Threshold) %>%
    summarise(n_regions = n_distinct(Region), .groups = "drop") %>%
    filter(n_regions == 2)
  
  if (nrow(metrics_both) > 0) {
    
    both_data <- all_results %>%
      semi_join(metrics_both, by = c("Metric_Label", "Lag_Label", 
                                     "HomeRangeSize", "Threshold"))
    
    both_forest <- ggplot(both_data,
                          aes(x = HomeRangeSize, y = Coefficient,
                              color = Region, shape = Region)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = Coefficient - StdError,
                        ymax = Coefficient + StdError),
                    position = position_dodge(width = 0.5),
                    width = 0.3,
                    linewidth = 0.7) +
      facet_wrap(~ Metric_Label + Lag_Label, scales = 'free', ncol = 4) +
      scale_color_manual(values = c("Bransfield" = "#2166AC",
                                    "Central_WAP" = "#B2182B")) +
      scale_shape_manual(values = c("Bransfield" = 16,
                                    "Central_WAP" = 17)) +
      labs(title = "Metrics Significant in BOTH Regions",
           subtitle = "Direct comparison of effect sizes where both regions show significance",
           x = "Home Range Size",
           y = "Effect Size (Coefficient ± SE)",
           color = "Region",
           shape = "Region") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
    
    ggsave(file.path(forest_dir, "Regional_Forest_Both_Regions.png"),
           both_forest, width = 14, height = 10, dpi = 300)
    
    ggsave(file.path(forest_dir, "Regional_Forest_Both_Regions.pdf"),
           both_forest, width = 14, height = 10)
    
    cat("  Saved: Regional_Forest_Both_Regions.png/pdf\n")
  } else {
    cat("  No metrics significant in both regions\n")
  }
  
  # ========================================================================
  # PLOT 4: Simplified forest plot by Lag only (averaged across home ranges)
  # ========================================================================
  cat("\nCreating simplified forest plot (averaged by lag)...\n")
  
  lag_summary <- all_results %>%
    group_by(Region, Metric_Label, Lag_Label) %>%
    summarise(
      Mean_Coefficient = mean(Coefficient, na.rm = TRUE),
      SE_Coefficient = sqrt(mean(StdError^2, na.rm = TRUE)),  # Pool SEs
      N_Models = n(),
      .groups = "drop"
    )
  
  lag_forest <- ggplot(lag_summary,
                       aes(x = Lag_Label, y = Mean_Coefficient,
                           color = Region, shape = Region)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_point(position = position_dodge(width = 0.5), size = 3.5) +
    geom_errorbar(aes(ymin = Mean_Coefficient - SE_Coefficient,
                      ymax = Mean_Coefficient + SE_Coefficient),
                  position = position_dodge(width = 0.5),
                  width = 0.3,
                  linewidth = 0.8) +
    facet_wrap(~ Metric_Label, scales = 'free_y', ncol = 4) +
    scale_color_manual(values = c("Bransfield" = "#2166AC",
                                  "Central_WAP" = "#B2182B")) +
    scale_shape_manual(values = c("Bransfield" = 16,
                                  "Central_WAP" = 17)) +
    labs(title = "Regional Comparison: Average Effect Sizes by Lag",
         subtitle = "Mean coefficients across all home ranges",
         x = "Lag Period",
         y = "Mean Effect Size (Coefficient ± SE)",
         color = "Region",
         shape = "Region") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(forest_dir, "Regional_Forest_Simplified_by_Lag.png"),
         lag_forest, width = 14, height = 8, dpi = 300)
  
  ggsave(file.path(forest_dir, "Regional_Forest_Simplified_by_Lag.pdf"),
         lag_forest, width = 14, height = 8)
  
  cat("  Saved: Regional_Forest_Simplified_by_Lag.png/pdf\n")
  
  cat("\n", rep("=", 70), "\n")
  cat("FOREST PLOT GENERATION COMPLETE\n")
  cat(rep("=", 70), "\n")
  
  return(all_results)
}