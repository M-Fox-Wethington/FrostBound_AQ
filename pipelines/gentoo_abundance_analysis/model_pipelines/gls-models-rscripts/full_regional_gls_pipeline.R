# ============================================================================
# GLS Analysis Pipeline - Regional and Pooled Options
# Gentoo Penguin Growth Rates vs Sea Ice Metrics
# ============================================================================

library(nlme)
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(showtext)
library(patchwork)

# Suppress showtext warnings
suppressMessages({
  font <- "Gudea"
  font_add_google(family = font, font, db_cache = TRUE)
  showtext_auto()
})

theme_set(theme_minimal(base_family = font, base_size = 12))

# ============================================================================
# SOURCE HELPER MODULES
# ============================================================================

analysis_root_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/gentoo_abundance_model_pipeline/gls-models-rscripts"

source(file.path(analysis_root_dir, "R/data_processing_helpers.R"))
source(file.path(analysis_root_dir, "R/model_helpers.R"))
source(file.path(analysis_root_dir, "R/visualization_helpers.R"))
source(file.path(analysis_root_dir, "R/interaction_helpers.R"))
source(file.path(analysis_root_dir, "R/pooling_helpers.R"))

# ============================================================================
# POOLED ANALYSIS FUNCTION (No Regional Split)
# ============================================================================

run_pooled_gls_analysis <- function(penguin_abundance_data,
                                    metric_files,
                                    results_dir,
                                    lags = 1:5,
                                    min_rows_for_fit = 30,
                                    metrics_priority = c(
                                      "mean_sic", "sd_sic", "ice_extent_km2",
                                      "mean_duration", "sd_duration",
                                      "mean_persistence", "sd_persistence"
                                    )) {
  
  cat("\n", rep("=", 70), "\n")
  cat("POOLED GLS ANALYSIS (All Regions Combined)\n")
  cat(rep("=", 70), "\n\n")
  
  # Create results directory
  pooled_dir <- file.path(results_dir, "Pooled_Analysis")
  dir.create(pooled_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Prepare penguin data (no regional filtering)
  penguin_data <- penguin_abundance_data %>%
    mutate(year = 1970 + season - 1) %>%
    filter(growth_rate <= 3, !is.na(growth_rate))
  
  cat(sprintf("Total colonies: %d\n", n_distinct(penguin_data$site_id)))
  cat(sprintf("Total observations: %d\n\n", nrow(penguin_data)))
  
  # Initialize results storage
  all_results_list <- list()
  
  # LOOP THROUGH EACH METRIC FILE
  for (f in metric_files) {
    cat(sprintf("\nProcessing file: %s\n", basename(f)))
    
    dat <- fread(f)
    
    # Normalize presence of grouping columns
    if (!"HomeRangeSize" %in% names(dat)) dat[, HomeRangeSize := "Unknown"]
    if (!"Threshold" %in% names(dat)) dat[, Threshold := NA_real_]
    
    # Which metric columns are in THIS file?
    metric_cols <- intersect(metrics_priority, names(dat))
    if (length(metric_cols) == 0) {
      warning(sprintf("  No recognized metric columns in file: %s", basename(f)))
      next
    }
    
    cat(sprintf("  Found metrics: %s\n", paste(metric_cols, collapse = ", ")))
    
    home_sizes <- sort(unique(dat$HomeRangeSize))
    thresholds <- sort(unique(dat$Threshold))
    
    # LOOP THROUGH HOME RANGE SIZES AND THRESHOLDS
    for (h in home_sizes) {
      for (th in thresholds) {
        
        df_ht <- dat %>% filter(HomeRangeSize == h, is.na(Threshold) | Threshold == th)
        if (nrow(df_ht) == 0) next
        
        cat(sprintf("  HomeRange: %s, Threshold: %s\n", h, as.character(th)))
        
        # LOOP THROUGH EACH METRIC COLUMN
        for (metric_col in metric_cols) {
          
          # Map metric_col to a friendly name
          metric_name <- dplyr::case_when(
            metric_col == "mean_sic"         ~ "SIC",
            metric_col == "sd_sic"           ~ "SIC_SD",
            metric_col == "ice_extent_km2"   ~ "Extent",
            metric_col == "mean_duration"    ~ "Duration",
            metric_col == "sd_duration"      ~ "Duration_SD",
            metric_col == "mean_persistence" ~ "Persistence",
            metric_col == "sd_persistence"   ~ "Persistence_SD",
            TRUE                              ~ metric_col
          )
          
          cat(sprintf("    Metric: %s (%s)\n", metric_name, metric_col))
          
          # Build yearly time series for this metric
          d_year <- tryCatch(
            yearly_aggregate_metric(df_ht, metric_col),
            error = function(e) {
              warning(sprintf("      Skipping %s: %s", metric_col, e$message))
              return(NULL)
            }
          )
          if (is.null(d_year)) next
          
          # Fit GLS models across lags (pooled data)
          gls_results <- fit_gls_models_pooled(penguin_data, d_year, metric_name)
          
          # Extract results to data frames
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
  
  # Export compiled results
  if (length(all_results_list) > 0) {
    all_results_df <- bind_rows(all_results_list)
    
    # Save all results
    write.csv(all_results_df,
              file.path(pooled_dir, "pooled_model_results_all_metrics.csv"),
              row.names = FALSE)
    cat(sprintf("\nSaved all results: %s\n",
                file.path(pooled_dir, "pooled_model_results_all_metrics.csv")))
    
    # Save significant results
    significant_results <- all_results_df %>%
      filter(Bonferroni_Adjusted_pValue < 0.05)
    
    if (nrow(significant_results) > 0) {
      write.csv(significant_results,
                file.path(pooled_dir, "pooled_model_results_significant.csv"),
                row.names = FALSE)
      cat(sprintf("Saved significant results: %s\n",
                  file.path(pooled_dir, "pooled_model_results_significant.csv")))
      
      # Summary statistics
      cat("\n=== POOLED ANALYSIS SUMMARY ===\n")
      cat(sprintf("Total significant results: %d\n", nrow(significant_results)))
      cat(sprintf("Unique metrics with significance: %d\n",
                  n_distinct(significant_results$Metric)))
      cat(sprintf("Mean N observations: %.0f\n", mean(significant_results$N_Obs)))
      cat(sprintf("Mean N colonies: %.0f\n", mean(significant_results$N_Colonies)))
    } else {
      cat("No significant results found\n")
    }
  }
  
  cat("\n", rep("=", 70), "\n")
  cat("POOLED ANALYSIS COMPLETE\n")
  cat(rep("=", 70), "\n")
  
  return(invisible(NULL))
}

# ============================================================================
# REGIONAL ANALYSIS FUNCTION (From existing pipeline)
# ============================================================================

run_regional_gls_analysis <- function(penguin_abundance_data, 
                                      metric_files,
                                      results_dir,
                                      central_cutoff_lat = -63.2,
                                      regions = c("Bransfield", "Central_WAP"),
                                      lags = 1:5,
                                      min_rows_for_fit = 30,
                                      metrics_priority = c(
                                        "mean_sic", "sd_sic", "ice_extent_km2",
                                        "mean_duration", "sd_duration",
                                        "mean_persistence", "sd_persistence"
                                      )) {
  
  # Create main results directory
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Identify latitude column
  lat_col <- pick_lat_col(penguin_abundance_data)
  if (is.null(lat_col)) {
    stop("No latitude column found in penguin data.")
  }
  
  # Add year and region columns to penguin data
  penguin_data_full <- penguin_abundance_data %>%
    mutate(
      year = 1970 + season - 1,
      .lat_raw = .data[[lat_col]],
      .lat_sgn = if_else(.lat_raw < 0, .lat_raw, -abs(.lat_raw)),
      region = if_else(.lat_sgn <= central_cutoff_lat, "Central_WAP", "Bransfield")
    ) %>%
    filter(growth_rate <= 3, !is.na(growth_rate))
  
  # LOOP THROUGH EACH REGION
  for (region_name in regions) {
    
    cat("\n", rep("=", 70), "\n")
    cat(sprintf("ANALYZING REGION: %s\n", region_name))
    cat(rep("=", 70), "\n\n")
    
    # Create region-specific directory
    region_dir <- file.path(results_dir, region_name)
    dir.create(region_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Subset penguin data for this region
    penguin_data_region <- penguin_data_full %>%
      filter(region == region_name)
    
    cat(sprintf("  N colonies in %s: %d\n", region_name, n_distinct(penguin_data_region$site_id)))
    cat(sprintf("  N observations in %s: %d\n\n", region_name, nrow(penguin_data_region)))
    
    # Initialize results storage
    all_results_list <- list()
    
    # LOOP THROUGH EACH METRIC FILE
    for (f in metric_files) {
      cat(sprintf("\nProcessing file: %s\n", basename(f)))
      
      dat <- fread(f)
      
      # Normalize presence of grouping columns
      if (!"HomeRangeSize" %in% names(dat)) dat[, HomeRangeSize := "Unknown"]
      if (!"Threshold" %in% names(dat)) dat[, Threshold := NA_real_]
      
      # Which metric columns are in THIS file?
      metric_cols <- intersect(metrics_priority, names(dat))
      if (length(metric_cols) == 0) {
        warning(sprintf("  No recognized metric columns in file: %s", basename(f)))
        next
      }
      
      cat(sprintf("  Found metrics: %s\n", paste(metric_cols, collapse = ", ")))
      
      home_sizes <- sort(unique(dat$HomeRangeSize))
      thresholds <- sort(unique(dat$Threshold))
      
      # LOOP THROUGH HOME RANGE SIZES AND THRESHOLDS
      for (h in home_sizes) {
        for (th in thresholds) {
          
          df_ht <- dat %>% filter(HomeRangeSize == h, is.na(Threshold) | Threshold == th)
          if (nrow(df_ht) == 0) next
          
          cat(sprintf("  HomeRange: %s, Threshold: %s\n", h, as.character(th)))
          
          # LOOP THROUGH EACH METRIC COLUMN
          for (metric_col in metric_cols) {
            
            # Map metric_col to a friendly name
            metric_name <- dplyr::case_when(
              metric_col == "mean_sic"         ~ "SIC",
              metric_col == "sd_sic"           ~ "SIC_SD",
              metric_col == "ice_extent_km2"   ~ "Extent",
              metric_col == "mean_duration"    ~ "Duration",
              metric_col == "sd_duration"      ~ "Duration_SD",
              metric_col == "mean_persistence" ~ "Persistence",
              metric_col == "sd_persistence"   ~ "Persistence_SD",
              TRUE                              ~ metric_col
            )
            
            cat(sprintf("    Metric: %s (%s)\n", metric_name, metric_col))
            
            # Build yearly time series for this metric
            d_year <- tryCatch(
              yearly_aggregate_metric(df_ht, metric_col),
              error = function(e) {
                warning(sprintf("      Skipping %s: %s", metric_col, e$message))
                return(NULL)
              }
            )
            if (is.null(d_year)) next
            
            # Fit GLS models across lags
            gls_results <- fit_gls_models_indiv(penguin_data_region, d_year, metric_name)
            
            # Extract results to data frames
            if (length(gls_results) > 0) {
              for (lag_name in names(gls_results)) {
                result <- gls_results[[lag_name]]
                
                result_df <- data.frame(
                  Region = region_name,
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
                  stringsAsFactors = FALSE
                )
                
                all_results_list[[length(all_results_list) + 1]] <- result_df
              }
            }
          }
        }
      }
    }
    
    # Export compiled results for this region
    if (length(all_results_list) > 0) {
      all_results_df <- bind_rows(all_results_list)
      
      # Save all results
      write.csv(all_results_df, 
                file.path(region_dir, "model_results_all_metrics.csv"), 
                row.names = FALSE)
      cat(sprintf("\nSaved all results: %s\n", 
                  file.path(region_dir, "model_results_all_metrics.csv")))
      
      # Save significant results
      significant_results <- all_results_df %>%
        filter(Bonferroni_Adjusted_pValue < 0.05)
      
      if (nrow(significant_results) > 0) {
        write.csv(significant_results, 
                  file.path(region_dir, "model_results_significant.csv"), 
                  row.names = FALSE)
        cat(sprintf("Saved significant results: %s\n", 
                    file.path(region_dir, "model_results_significant.csv")))
      } else {
        cat(sprintf("No significant results for %s\n", region_name))
      }
    }
  }
  
  cat("\n", rep("=", 70), "\n")
  cat("REGIONAL ANALYSIS COMPLETE\n")
  cat(rep("=", 70), "\n")
}

# ============================================================================
# MAIN PIPELINE EXECUTION
# ============================================================================

# Define file paths
penguin_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/data/gentoo-abundance-model/inputs/modeled_gentoo_parameters.csv"

metric_files <- c(
  "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/data/gentoo-abundance-model/metric-calculation-csv/daily_sic_statistics.csv",
  "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/data/gentoo-abundance-model/metric-calculation-csv/sea_ice_duration_persistence_stats.csv")

results_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/results"

# Load penguin data
penguin_data <- fread(penguin_path)

# ============================================================================
# OPTION 1: RUN REGIONAL ANALYSIS FIRST (for comparison/justification)
# ============================================================================

cat("\n### STEP 1: REGIONAL ANALYSIS ###\n")
run_regional_gls_analysis(
  penguin_abundance_data = penguin_data,
  metric_files = metric_files,
  results_dir = results_path,
  central_cutoff_lat = -63.2,
  regions = c("Bransfield", "Central_WAP"),
  lags = 1:5,
  min_rows_for_fit = 30
)

# ============================================================================
# OPTION 2: TEST REGIONAL INTERACTIONS
# ============================================================================

cat("\n### STEP 2: INTERACTION ANALYSIS ###\n")
interaction_results <- test_regional_interactions(
  penguin_abundance_data = penguin_data,
  metric_files = metric_files,
  results_dir = results_path,
  central_cutoff_lat = -63.2,
  lags = 1:5,
  min_rows_for_fit = 30
)

# ============================================================================
# OPTION 3: ANALYZE POOLING JUSTIFICATION
# ============================================================================

cat("\n### STEP 3: POOLING JUSTIFICATION ANALYSIS ###\n")
pooling_analysis <- analyze_pooling_justification(
  regional_comparison_path = file.path(results_path, "Regional_Comparison_Summary.csv"),
  interaction_results_path = file.path(results_path, "Interaction_Analysis", "interaction_tests_all.csv"),
  output_dir = file.path(results_path, "Pooling_Justification")
)

# ============================================================================
# OPTION 4: RUN POOLED ANALYSIS (Primary Analysis)
# ============================================================================

if (pooling_analysis$recommendation == "POOL REGIONS") {
  cat("\n### STEP 4: POOLED ANALYSIS (PRIMARY) ###\n")
  run_pooled_gls_analysis(
    penguin_abundance_data = penguin_data,
    metric_files = metric_files,
    results_dir = results_path,
    lags = 1:5,
    min_rows_for_fit = 30
  )
} else {
  cat("\n### POOLING NOT RECOMMENDED - Use regional results as primary ###\n")
}

# ============================================================================
# OPTION 5: VISUALIZATION (Regional comparison for supplementary)
# ============================================================================

cat("\n### STEP 5: VISUALIZATION ###\n")

# Regional comparison plots (for supplement)
tryCatch({
  create_regional_comparison_plots(
    results_dir = results_path,
    penguin_data_path = penguin_path,
    metric_files = metric_files,
    central_cutoff_lat = -63.2
  )
}, error = function(e) {
  cat("\nError during plot generation:", e$message, "\n")
})

# Forest plots (regional comparison)
forest_results <- create_regional_forest_plots(results_path)

# Summary table
regional_comparison <- create_regional_comparison_table(results_path)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n\n")
cat(rep("=", 70), "\n")
cat("GLS ANALYSIS PIPELINE COMPLETE!\n")
cat(rep("=", 70), "\n")
cat("\nOutputs saved in:", results_path, "\n")
cat("  - Pooled_Analysis/: PRIMARY RESULTS (pooled across regions)\n")
cat("  - Bransfield/: Regional results (supplementary)\n")
cat("  - Central_WAP/: Regional results (supplementary)\n")
cat("  - Pooling_Justification/: Statistical support for pooling\n")
cat("  - Interaction_Analysis/: Formal interaction tests\n")
cat("  - Comparison_Plots/: Side-by-side visualizations\n")
cat("  - Forest_Plots/: Regional comparison forest plots\n")
cat(sprintf("\nFinal Recommendation: %s\n", pooling_analysis$recommendation))

