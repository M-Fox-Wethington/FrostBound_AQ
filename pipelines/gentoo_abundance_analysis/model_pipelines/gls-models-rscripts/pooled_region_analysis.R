# ============================================================================
# COMPREHENSIVE POOLED ANALYSIS WITH VISUALIZATIONS
# Gentoo Penguin Growth Rates vs Sea Ice Metrics
# All Regions Combined (Bransfield + Central WAP)
# ============================================================================

# ── LIBRARIES ────────────────────────────────────────────────────────────────
library(nlme)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(viridis)
library(patchwork)

# ── CONFIGURATION ────────────────────────────────────────────────────────────
# Define file paths (UPDATE THESE TO YOUR PATHS)
penguin_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/data/inputs/modeled_gentoo_parameters.csv"
metric_files <- c(
  "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/data/gentoo-abundance-model/metric-calculation-csv/daily_sic_statistics.csv",
  "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/data/gentoo-abundance-model/metric-calculation-csv/sea_ice_duration_persistence_stats.csv"
)

results_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/results/Pooled_Analysis"

# Metrics to analyze
metrics_priority <- c(
  "mean_sic", "sd_sic", "ice_extent_km2",
  "mean_duration", "sd_duration",
  "mean_persistence", "sd_persistence"
)

# ============================================================================
# HELPER FUNCTION: Get individual ice metrics
# ============================================================================
get_individual_ice_metrics <- function(year, ice_data, metric) {
  ice_values <- ice_data %>%
    filter(year == !!year) %>%
    pull(!!sym(metric))
  
  if (length(ice_values) > 0) {
    return(mean(ice_values, na.rm = TRUE))
  } else {
    return(NA)
  }
}

# ============================================================================
# HELPER FUNCTION: Yearly aggregate metric
# ============================================================================
yearly_aggregate_metric <- function(dat, metric_col) {
  if (!metric_col %in% names(dat)) {
    stop(sprintf("Metric column '%s' not found", metric_col))
  }
  
  # Convert date to year if date column exists
  if ("date" %in% names(dat) && !"year" %in% names(dat)) {
    dat <- dat %>%
      mutate(year = year(as.Date(date)))
  }
  
  # Check if year column exists now
  if (!"year" %in% names(dat)) {
    stop("No 'year' or 'date' column found in data")
  }
  
  d_year <- dat %>%
    group_by(year) %>%
    summarise(
      metric_value = mean(get(metric_col), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(metric_value))
  
  if (nrow(d_year) == 0) {
    stop("No valid yearly data after aggregation")
  }
  
  return(d_year)
}


# ============================================================================
# HELPER FUNCTION: Order home range sizes numerically
# ============================================================================
order_home_range <- function(hr_vector) {
  # Extract numeric values from strings like "25km", "100km", etc.
  numeric_values <- as.numeric(gsub("km", "", hr_vector))
  
  # Create factor with levels ordered by numeric value
  factor(hr_vector, levels = unique(hr_vector[order(numeric_values)]))
}

# ============================================================================
# CORE FUNCTION: Fit GLS models for pooled data
# ============================================================================
fit_gls_models_pooled <- function(penguin_data, ice_data, metric) {
  results <- list()
  
  for (lag in 1:5) {
    penguin_lagged <- penguin_data %>%
      rowwise() %>%
      mutate(overwinter_ice = get_individual_ice_metrics(year - lag, ice_data, metric)) %>%
      ungroup() %>%
      filter(!is.na(overwinter_ice))
    
    if (nrow(penguin_lagged) < 30) {
      cat(sprintf("    Skipping lag %d: insufficient data (n=%d)\n", lag, nrow(penguin_lagged)))
      next
    }
    
    cat(sprintf("    Fitting GLS model for lag %d (n=%d)\n", lag, nrow(penguin_lagged)))
    
    tryCatch({
      formula <- as.formula("growth_rate ~ overwinter_ice")
      model <- gls(formula, 
                   correlation = corAR1(form = ~1 | site_id), 
                   data = penguin_lagged)
      
      summary_model <- summary(model)
      bonferroni_p_value <- p.adjust(summary_model$tTable[2, 4], method = "bonferroni", n = 5)
      
      results[[paste("Lag", lag)]] <- list(
        lag = lag,
        AIC = AIC(model),
        BIC = BIC(model),
        coefficients = summary_model$tTable,
        p_value = summary_model$tTable[2, 4],
        bonferroni_p_value = bonferroni_p_value,
        is_significant = bonferroni_p_value < 0.05,
        model = model,
        penguin_data = penguin_lagged
      )
    }, error = function(e) {
      cat(sprintf("    Error fitting model for lag %d: %s\n", lag, e$message))
    })
  }
  
  return(results)
}

# ============================================================================
# MAIN FUNCTION: Run pooled GLS analysis
# ============================================================================
run_pooled_analysis <- function(penguin_data, 
                                metric_files,
                                results_dir,
                                lags = 1:5,
                                min_rows_for_fit = 30) {
  
  cat("\n", rep("=", 70), "\n")
  cat("POOLED GLS ANALYSIS (ALL REGIONS COMBINED)\n")
  cat(rep("=", 70), "\n\n")
  
  # Create results directory
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Prepare penguin data (no regional filtering)
  penguin_full <- penguin_data %>%
    mutate(year = 1970 + season - 1) %>%
    filter(growth_rate <= 3, !is.na(growth_rate))
  
  cat(sprintf("Total colonies: %d\n", n_distinct(penguin_full$site_id)))
  cat(sprintf("Total observations: %d\n\n", nrow(penguin_full)))
  
  all_results_list <- list()
  
  # Loop through each metric file
  for (f in metric_files) {
    cat(sprintf("\nProcessing file: %s\n", basename(f)))
    
    dat <- fread(f)
    
    # Normalize columns
    if (!"HomeRangeSize" %in% names(dat)) dat[, HomeRangeSize := "Unknown"]
    if (!"Threshold" %in% names(dat)) dat[, Threshold := NA_real_]
    
    # Identify metric columns
    metric_cols <- intersect(metrics_priority, names(dat))
    if (length(metric_cols) == 0) {
      cat(sprintf("  No recognized metrics in %s\n", basename(f)))
      next
    }
    
    cat(sprintf("  Found metrics: %s\n", paste(metric_cols, collapse = ", ")))
    
    home_sizes <- sort(unique(dat$HomeRangeSize))
    thresholds <- sort(unique(dat$Threshold))
    
    # Loop through combinations
    for (h in home_sizes) {
      for (th in thresholds) {
        
        df_ht <- dat %>% 
          filter(HomeRangeSize == h, is.na(Threshold) | Threshold == th)
        
        if (nrow(df_ht) == 0) next
        
        cat(sprintf("\n  HomeRange: %s, Threshold: %s\n", h, as.character(th)))
        
        # Loop through each metric
        for (metric_col in metric_cols) {
          
          # Map to friendly name
          metric_name <- case_when(
            metric_col == "mean_sic"         ~ "SIC",
            metric_col == "sd_sic"           ~ "SIC Variability",
            metric_col == "ice_extent_km2"   ~ "Sea Ice Extent",
            metric_col == "mean_duration"    ~ "Duration",
            metric_col == "sd_duration"      ~ "Duration Variability",
            metric_col == "mean_persistence" ~ "Open Water Frequency",
            metric_col == "sd_persistence"   ~ "Open Water Frequency Variability",
            TRUE                              ~ metric_col
          )
          
          cat(sprintf("    Metric: %s (%s)\n", metric_name, metric_col))
          
          # Build yearly time series
          d_year <- tryCatch(
            yearly_aggregate_metric(df_ht, metric_col),
            error = function(e) {
              cat(sprintf("      Error: %s\n", e$message))
              return(NULL)
            }
          )
          if (is.null(d_year)) next
          
          # Fit GLS models
          gls_results <- fit_gls_models_pooled(penguin_full, d_year, "metric_value")
          
          # Extract results
          if (length(gls_results) > 0) {
            for (lag_name in names(gls_results)) {
              result <- gls_results[[lag_name]]
              
              result_df <- data.frame(
                Analysis = "Pooled",
                File = basename(f),
                Metric = metric_name,
                MetricColumn = metric_col,
                HomeRangeSize = h,
                Threshold = th,
                Lag = result$lag,
                N_Obs = nrow(result$penguin_data),
                N_Colonies = n_distinct(result$penguin_data$site_id),
                AIC = result$AIC,
                BIC = result$BIC,
                Coefficient_Intercept = result$coefficients[1, 1],
                StdError_Intercept = result$coefficients[1, 2],
                tValue_Intercept = result$coefficients[1, 3],
                pValue_Intercept = result$coefficients[1, 4],
                Coefficient = result$coefficients[2, 1],
                StdError = result$coefficients[2, 2],
                tValue = result$coefficients[2, 3],
                pValue = result$coefficients[2, 4],
                Bonferroni_Adjusted_pValue = result$bonferroni_p_value,
                Is_Significant = result$is_significant,
                stringsAsFactors = FALSE
              )
              
              all_results_list[[length(all_results_list) + 1]] <- result_df
            }
          }
        }
      }
    }
  }
  
  # Save results
  if (length(all_results_list) > 0) {
    all_results_df <- bind_rows(all_results_list)
    
    # Save all results
    write.csv(all_results_df,
              file.path(results_dir, "pooled_model_results_all_metrics.csv"),
              row.names = FALSE)
    cat(sprintf("\nSaved: %s\n", file.path(results_dir, "pooled_model_results_all_metrics.csv")))
    
    # Save significant results
    significant_results <- all_results_df %>%
      filter(Is_Significant == TRUE)
    
    if (nrow(significant_results) > 0) {
      write.csv(significant_results,
                file.path(results_dir, "pooled_model_results_significant.csv"),
                row.names = FALSE)
      cat(sprintf("Saved: %s\n", file.path(results_dir, "pooled_model_results_significant.csv")))
      
      cat("\n=== POOLED ANALYSIS SUMMARY ===\n")
      cat(sprintf("Total significant results: %d\n", nrow(significant_results)))
      cat(sprintf("Unique metrics: %d\n", n_distinct(significant_results$Metric)))
      cat(sprintf("Home range sizes: %d\n", n_distinct(significant_results$HomeRangeSize)))
      cat(sprintf("Lag periods: %d\n\n", n_distinct(significant_results$Lag)))
      
      # Metric-level summary
      metric_summary <- significant_results %>%
        group_by(Metric) %>%
        summarise(
          N_Significant = n(),
          Mean_Coefficient = mean(Coefficient),
          Range_Coefficient = paste0("[", 
                                     round(min(Coefficient), 3), ", ",
                                     round(max(Coefficient), 3), "]"),
          .groups = "drop"
        )
      
      print(metric_summary)
    } else {
      cat("No significant results found.\n")
    }
    
    return(all_results_df)
  } else {
    cat("No results generated.\n")
    return(NULL)
  }
}

# ============================================================================
# VISUALIZATION 1: Forest Plots - Pooled Results
# ============================================================================
create_pooled_forest_plots <- function(results_dir) {
  
  cat("\n", rep("=", 70), "\n")
  cat("CREATING POOLED FOREST PLOTS\n")
  cat(rep("=", 70), "\n\n")
  
  forest_dir <- file.path(results_dir, "Forest_Plots")
  dir.create(forest_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load significant results
  sig_path <- file.path(results_dir, "pooled_model_results_significant.csv")
  
  if (!file.exists(sig_path)) {
    cat("No significant results found at:", sig_path, "\n")
    return(NULL)
  }
  
  dat <- read.csv(sig_path, stringsAsFactors = FALSE)
  
  if (nrow(dat) == 0) {
    cat("Significant results file is empty.\n")
    return(NULL)
  }
  
  cat(sprintf("Loaded %d significant results\n", nrow(dat)))
  
  # Prepare data for plotting
  plot_data <- dat %>%
    mutate(
      hr_size = order_home_range(HomeRangeSize),
      lag = factor(Lag, levels = as.character(1:5)),
      threshold = if_else(is.na(Threshold), "NA", as.character(Threshold)),
      metric = Metric,
      estimate = Coefficient
    ) %>%
    # Separate main vs variability metrics
    mutate(
      metric_type = if_else(str_detect(metric, "Variability"), "Variability", "Main")
    )
  
  # Complete grid for missing combinations
  plot_data <- plot_data %>%
    complete(metric, threshold, hr_size, lag)
  
  # Separate datasets
  main_df <- plot_data %>% filter(metric_type == "Main" | is.na(metric_type))
  var_df <- plot_data %>% filter(metric_type == "Variability" | is.na(metric_type))
  
  # ── PLOT 1: Main metrics forest plot ────────────────────────────────────
  cat("\nCreating main metrics forest plot...\n")
  
  if (nrow(main_df %>% filter(!is.na(estimate))) > 0) {
    
    forest_main <- main_df %>%
      filter(!is.na(estimate)) %>%
      ggplot(aes(x = estimate, y = hr_size, color = lag)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(
        aes(xmin = estimate - 1.96 * StdError,
            xmax = estimate + 1.96 * StdError),
        height = 0.2, 
        position = position_dodge(width = 0.7),
        na.rm = TRUE
      ) +
      geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
      facet_grid(metric ~ threshold, scales = "free_y") +
      scale_color_viridis_d(name = "Lag (years)", end = 0.8) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      ) +
      labs(
        title = "Pooled Analysis: Main Sea Ice Metrics",
        x = "Coefficient (95% CI)",
        y = "Home Range Size (km)"
      )
    
    print(forest_main)
    ggsave(file.path(forest_dir, "forest_main_metrics.png"), 
           forest_main, width = 10, height = 8, dpi = 300)
    ggsave(file.path(forest_dir, "forest_main_metrics.pdf"), 
           forest_main, width = 10, height = 8)
    
    cat("  Saved: forest_main_metrics.png/pdf\n")
  }
  
  # ── PLOT 2: Variability metrics forest plot ─────────────────────────────
  cat("\nCreating variability metrics forest plot...\n")
  
  if (nrow(var_df %>% filter(!is.na(estimate))) > 0) {
    
    forest_var <- var_df %>%
      filter(!is.na(estimate)) %>%
      ggplot(aes(x = estimate, y = hr_size, color = lag)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(
        aes(xmin = estimate - 1.96 * StdError,
            xmax = estimate + 1.96 * StdError),
        height = 0.2, 
        position = position_dodge(width = 0.7),
        na.rm = TRUE
      ) +
      geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
      facet_grid(metric ~ threshold, scales = "free_y") +
      scale_color_viridis_d(name = "Lag (years)", end = 0.8) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      ) +
      labs(
        title = "Pooled Analysis: Variability Metrics",
        x = "Coefficient (95% CI)",
        y = "Home Range Size (km)"
      )
    
    print(forest_var)
    ggsave(file.path(forest_dir, "forest_variability_metrics.png"), 
           forest_var, width = 10, height = 8, dpi = 300)
    ggsave(file.path(forest_dir, "forest_variability_metrics.pdf"), 
           forest_var, width = 10, height = 8)
    
    cat("  Saved: forest_variability_metrics.png/pdf\n")
  }
  
  # ── PLOT 3: Per-metric forest plots ─────────────────────────────────────
  cat("\nCreating individual metric forest plots...\n")
  
  # Derive base metric names
  base_metrics <- plot_data %>%
    pull(metric) %>%
    str_remove(" Variability$") %>%
    unique() %>%
    na.omit()
  
  for (bm in base_metrics) {
    df_sub <- plot_data %>%
      filter(metric %in% c(bm, paste0(bm, " Variability"))) %>%
      filter(!is.na(estimate))
    
    if (nrow(df_sub) == 0) next
    
    p <- ggplot(df_sub, aes(x = estimate, y = hr_size, color = lag)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(
        aes(xmin = estimate - 1.96 * StdError,
            xmax = estimate + 1.96 * StdError),
        height = 0.2, 
        position = position_dodge(width = 0.7),
        na.rm = TRUE
      ) +
      geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
      facet_grid(metric ~ threshold, scales = "free_y", switch = "y") +
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
        title = paste0("Forest Plot: ", bm, " & Variability"),
        x = "Coefficient (95% CI)",
        y = "Home Range Size (km)"
      )
    
    print(p)
    
    safe_name <- gsub("[^[:alnum:]_]", "_", bm)
    ggsave(
      filename = file.path(forest_dir, paste0("forest_", safe_name, ".png")),
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )
    ggsave(
      filename = file.path(forest_dir, paste0("forest_", safe_name, ".pdf")),
      plot = p,
      width = 10,
      height = 6
    )
    
    cat(sprintf("  Saved: forest_%s.png/pdf\n", safe_name))
  }
  
  cat("\n", rep("=", 70), "\n")
  cat("FOREST PLOT GENERATION COMPLETE\n")
  cat(rep("=", 70), "\n")
  
  return(plot_data)
}

# ============================================================================
# VISUALIZATION 2: Heatmaps - Pooled Results
# ============================================================================
create_pooled_heatmaps <- function(results_dir) {
  
  cat("\n", rep("=", 70), "\n")
  cat("CREATING POOLED HEATMAPS\n")
  cat(rep("=", 70), "\n\n")
  
  heatmap_dir <- file.path(results_dir, "Heatmaps")
  dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load significant results
  sig_path <- file.path(results_dir, "pooled_model_results_significant.csv")
  
  if (!file.exists(sig_path)) {
    cat("No significant results found.\n")
    return(NULL)
  }
  
  dat <- read.csv(sig_path, stringsAsFactors = FALSE)
  
  if (nrow(dat) == 0) {
    cat("Significant results file is empty.\n")
    return(NULL)
  }
  
  # Prepare data
  df2 <- dat %>%
    mutate(
      hr_size = order_home_range(HomeRangeSize),
      lag = factor(Lag, levels = as.character(1:5)),
      threshold = factor(Threshold, levels = sort(unique(Threshold), na.last = TRUE)),
      metric = Metric,
      estimate = Coefficient
    ) %>%
    complete(metric, threshold, hr_size, lag)
  
  # Separate main and variability
  main_df <- df2 %>%
    filter(!str_detect(metric, "Variability"))
  
  var_df <- df2 %>%
    filter(str_detect(metric, "Variability"))
  
  # ── HEATMAP 1: Main metrics ──────────────────────────────────────────────
  cat("\nCreating main metrics heatmap...\n")
  
  if (nrow(main_df %>% filter(!is.na(estimate))) > 0) {
    
    main_plot <- ggplot(main_df, aes(x = lag, y = hr_size, fill = estimate)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(
        aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
        size = 3.5, color = "white"
      ) +
      scale_fill_viridis(
        option = "D", 
        direction = -1,
        begin = 0.15, 
        end = 0.75,
        na.value = "grey95", 
        name = "Coef."
      ) +
      facet_grid(metric ~ threshold, scales = "free_y") +
      theme_minimal(base_size = 14, base_family = "sans") +
      theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "#F0F0F0", color = NA),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        legend.position = "right",
        legend.background = element_blank(),
        plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
      ) +
      labs(
        title = "Pooled Analysis: Main Sea Ice Metrics Coefficients",
        x = "Lag (years)",
        y = "Home Range Size (km)"
      )
    
    print(main_plot)
    ggsave(file.path(heatmap_dir, "heatmap_main_metrics.png"), 
           main_plot, width = 12, height = 8, dpi = 300)
    ggsave(file.path(heatmap_dir, "heatmap_main_metrics.pdf"), 
           main_plot, width = 12, height = 8)
    
    cat("  Saved: heatmap_main_metrics.png/pdf\n")
  }
  
  # ── HEATMAP 2: Variability metrics ───────────────────────────────────────
  cat("\nCreating variability metrics heatmap...\n")
  
  if (nrow(var_df %>% filter(!is.na(estimate))) > 0) {
    
    var_plot <- ggplot(var_df, aes(x = lag, y = hr_size, fill = estimate)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(
        aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
        size = 3.5, color = "white"
      ) +
      scale_fill_viridis(
        option = "D", 
        direction = -1,
        begin = 0.15, 
        end = 0.75,
        na.value = "grey95", 
        name = "Coef."
      ) +
      facet_grid(metric ~ threshold, scales = "free_y") +
      theme_minimal(base_size = 14, base_family = "sans") +
      theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "#F0F0F0", color = NA),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        legend.position = "right",
        legend.background = element_blank(),
        plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
      ) +
      labs(
        title = "Pooled Analysis: Variability Metrics Coefficients",
        x = "Lag (years)",
        y = "Home Range Size (km)"
      )
    
    print(var_plot)
    ggsave(file.path(heatmap_dir, "heatmap_variability_metrics.png"), 
           var_plot, width = 12, height = 8, dpi = 300)
    ggsave(file.path(heatmap_dir, "heatmap_variability_metrics.pdf"), 
           var_plot, width = 12, height = 8)
    
    cat("  Saved: heatmap_variability_metrics.png/pdf\n")
  }
  
  # ── HEATMAP 3: Combined full grid ────────────────────────────────────────
  cat("\nCreating combined full grid heatmap...\n")
  
  heat_all <- ggplot(df2, aes(x = lag, y = hr_size, fill = estimate)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
      size = 3.5, color = "white"
    ) +
    scale_fill_viridis(
      option = "D",
      direction = -1,
      begin = 0.15,
      end = 0.75,
      na.value = "grey95",
      name = "Coef."
    ) +
    facet_grid(metric ~ threshold, scales = "free_y") +
    theme_minimal(base_size = 14, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "#F0F0F0", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 13),
      legend.position = "right",
      legend.background = element_blank(),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
    ) +
    labs(
      title = "Effect Sizes of Winter Sea Ice Metrics on Gentoo Penguin Growth",
      subtitle = "Pooled Analysis: All Regions Combined",
      x = "Lag (years)",
      y = "Home Range Size (km)"
    )
  
  print(heat_all)
  ggsave(file.path(heatmap_dir, "heatmap_full_grid.png"),
         heat_all, width = 12, height = 8, dpi = 300)
  ggsave(file.path(heatmap_dir, "heatmap_full_grid.pdf"),
         heat_all, width = 12, height = 8)
  
  cat("  Saved: heatmap_full_grid.png/pdf\n")
  
  cat("\n", rep("=", 70), "\n")
  cat("HEATMAP GENERATION COMPLETE\n")
  cat(rep("=", 70), "\n")
  
  return(df2)
}

# ============================================================================
# VISUALIZATION 3: Summary Tables and Statistics
# ============================================================================
create_pooled_summary <- function(results_dir) {
  
  cat("\n", rep("=", 70), "\n")
  cat("CREATING POOLED SUMMARY STATISTICS\n")
  cat(rep("=", 70), "\n\n")
  
  # Load all results
  all_path <- file.path(results_dir, "pooled_model_results_all_metrics.csv")
  sig_path <- file.path(results_dir, "pooled_model_results_significant.csv")
  
  if (!file.exists(all_path)) {
    cat("No results found.\n")
    return(NULL)
  }
  
  all_results <- read.csv(all_path, stringsAsFactors = FALSE)
  
  # Summary by metric
  metric_summary <- all_results %>%
    group_by(Metric) %>%
    summarise(
      N_Models = n(),
      N_Significant = sum(Is_Significant),
      Percent_Significant = round(100 * N_Significant / N_Models, 1),
      Mean_Coefficient = mean(Coefficient, na.rm = TRUE),
      SD_Coefficient = sd(Coefficient, na.rm = TRUE),
      Range_Coefficient = paste0("[", 
                                 round(min(Coefficient, na.rm = TRUE), 4), ", ",
                                 round(max(Coefficient, na.rm = TRUE), 4), "]"),
      Mean_pValue = mean(pValue, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(N_Significant))
  
  write.csv(metric_summary,
            file.path(results_dir, "pooled_summary_by_metric.csv"),
            row.names = FALSE)
  
  cat("Saved: pooled_summary_by_metric.csv\n\n")
  print(metric_summary)
  
  # Summary by lag
  lag_summary <- all_results %>%
    group_by(Lag) %>%
    summarise(
      N_Models = n(),
      N_Significant = sum(Is_Significant),
      Percent_Significant = round(100 * N_Significant / N_Models, 1),
      Mean_Coefficient = mean(Coefficient, na.rm = TRUE),
      SD_Coefficient = sd(Coefficient, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Lag)
  
  write.csv(lag_summary,
            file.path(results_dir, "pooled_summary_by_lag.csv"),
            row.names = FALSE)
  
  cat("\nSaved: pooled_summary_by_lag.csv\n\n")
  print(lag_summary)
  
  # Summary by home range size
  hr_summary <- all_results %>%
    group_by(HomeRangeSize) %>%
    summarise(
      N_Models = n(),
      N_Significant = sum(Is_Significant),
      Percent_Significant = round(100 * N_Significant / N_Models, 1),
      Mean_Coefficient = mean(Coefficient, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(HomeRangeSize)
  
  write.csv(hr_summary,
            file.path(results_dir, "pooled_summary_by_homerange.csv"),
            row.names = FALSE)
  
  cat("\nSaved: pooled_summary_by_homerange.csv\n\n")
  print(hr_summary)
  
  cat("\n", rep("=", 70), "\n")
  
  return(list(
    metric_summary = metric_summary,
    lag_summary = lag_summary,
    hr_summary = hr_summary
  ))
}

# ============================================================================
# MAIN EXECUTION WORKFLOW
# ============================================================================

cat("\n")
cat(rep("=", 70), "\n")
cat("COMPREHENSIVE POOLED ANALYSIS PIPELINE\n")
cat("Gentoo Penguin Growth Rates vs Sea Ice Metrics\n")
cat("All Regions Combined (Bransfield + Central WAP)\n")
cat(rep("=", 70), "\n")

# Create main results directory
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

# Load penguin data
cat("\nLoading penguin abundance data...\n")
penguin_data <- fread(file = penguin_path)
cat(sprintf("Loaded %d observations from %d colonies\n", 
            nrow(penguin_data), n_distinct(penguin_data$site_id)))

# ── STEP 1: Pooled GLS Analysis ──────────────────────────────────────────────
cat("\n", rep("=", 70), "\n")
cat("STEP 1: POOLED GLS ANALYSIS\n")
cat(rep("=", 70), "\n")

pooled_results <- run_pooled_analysis(
  penguin_data = penguin_data,
  metric_files = metric_files,
  results_dir = results_path,
  lags = 1:5,
  min_rows_for_fit = 30
)

# ── STEP 2: Forest Plots ──────────────────────────────────────────────────────
cat("\n", rep("=", 70), "\n")
cat("STEP 2: FOREST PLOT GENERATION\n")
cat(rep("=", 70), "\n")

forest_plots <- tryCatch({
  create_pooled_forest_plots(results_path)
}, error = function(e) {
  cat("\nError creating forest plots:", e$message, "\n")
  return(NULL)
})

# ── STEP 3: Heatmaps ──────────────────────────────────────────────────────────
cat("\n", rep("=", 70), "\n")
cat("STEP 3: HEATMAP GENERATION\n")
cat(rep("=", 70), "\n")

heatmaps <- tryCatch({
  create_pooled_heatmaps(results_path)
}, error = function(e) {
  cat("\nError creating heatmaps:", e$message, "\n")
  return(NULL)
})

# ── STEP 4: Summary Statistics ────────────────────────────────────────────────
cat("\n", rep("=", 70), "\n")
cat("STEP 4: SUMMARY STATISTICS\n")
cat(rep("=", 70), "\n")

summary_stats <- tryCatch({
  create_pooled_summary(results_path)
}, error = function(e) {
  cat("\nError creating summary:", e$message, "\n")
  return(NULL)
})

# ── FINAL SUMMARY ─────────────────────────────────────────────────────────────
cat("\n\n")
cat(rep("=", 70), "\n")
cat("POOLED ANALYSIS PIPELINE COMPLETE!\n")
cat(rep("=", 70), "\n")
cat("\nOutputs saved in:", results_path, "\n")
cat("  - pooled_model_results_all_metrics.csv: All model results\n")
cat("  - pooled_model_results_significant.csv: Significant results only\n")
cat("  - Forest_Plots/: Forest plot visualizations\n")
cat("  - Heatmaps/: Heatmap visualizations\n")
cat("  - pooled_summary_*.csv: Summary statistics by metric/lag/homerange\n")
cat("\n", rep("=", 70), "\n")