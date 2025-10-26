# ============================================================================
# interaction_analysis_helpers.R
# Helper functions for testing regional differences in sea ice-penguin relationships
# 
# Dependencies: Requires the following functions to be loaded in the calling script:
#   - pick_lat_col()
#   - safe_year_from_date() 
#   - yearly_aggregate_metric()
#
# Required packages: nlme, data.table, dplyr, lubridate
# ============================================================================

#' Test for Significant Regional Differences Using Interaction Models
#'
#' This function tests whether sea ice effects on penguin growth differ 
#' significantly between Bransfield and Central WAP regions using 
#' metric × region interaction terms in GLS models.
#'
#' @param penguin_abundance_data Data frame with penguin growth rates
#' @param metric_files Character vector of paths to sea ice metric files
#' @param results_dir Directory path for saving results
#' @param central_cutoff_lat Latitude threshold for regional classification (default -63.2)
#' @param lags Integer vector of lag years to test (default 1:5)
#' @param min_rows_for_fit Minimum observations required for model fitting (default 30)
#'
#' @return Data frame with interaction test results, or NULL if no models fitted
#' @export
test_regional_interactions <- function(penguin_abundance_data,
                                       metric_files,
                                       results_dir,
                                       central_cutoff_lat = -63.2,
                                       lags = 1:5,
                                       min_rows_for_fit = 30) {
  
  cat("\n", rep("=", 70), "\n")
  cat("TESTING FOR SIGNIFICANT REGIONAL DIFFERENCES\n")
  cat(rep("=", 70), "\n\n")
  
  # Setup directories
  interaction_dir <- file.path(results_dir, "Interaction_Analysis")
  dir.create(interaction_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load and prepare penguin data with regions
  lat_col <- pick_lat_col(penguin_abundance_data)
  
  penguin_full <- penguin_abundance_data %>%
    mutate(
      year = 1970 + season - 1,
      .lat_raw = .data[[lat_col]],
      .lat_sgn = if_else(.lat_raw < 0, .lat_raw, -abs(.lat_raw)),
      region = if_else(.lat_sgn <= central_cutoff_lat, "Central_WAP", "Bransfield")
    ) %>%
    filter(growth_rate <= 3, !is.na(growth_rate))
  
  # Store all interaction test results
  interaction_results_list <- list()
  
  # Process each metric file
  for (f in metric_files) {
    cat(sprintf("\n--- Processing file: %s ---\n", basename(f)))
    
    # Load metric data
    metric_dat <- fread(f)
    
    # Identify which metric columns are present in this file
    potential_metric_cols <- c("mean_sic", "sd_sic", "ice_extent_km2", 
                               "mean_duration", "sd_duration",
                               "mean_persistence", "sd_persistence")
    
    metric_cols <- potential_metric_cols[potential_metric_cols %in% names(metric_dat)]
    
    if (length(metric_cols) == 0) {
      cat("  No recognized metric columns found. Skipping.\n")
      next
    }
    
    cat(sprintf("  Found %d metric columns: %s\n", 
                length(metric_cols), paste(metric_cols, collapse=", ")))
    
    # Get unique parameter combinations
    h_vals <- unique(metric_dat$HomeRangeSize)
    th_vals <- unique(metric_dat$Threshold)
    th_vals <- th_vals[!is.na(th_vals)]
    
    # Process each metric
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
      
      cat(sprintf("\n  Metric: %s (%s)\n", metric_name, metric_col))
      
      # Process each parameter combination
      for (h in h_vals) {
        for (th in th_vals) {
          
          # Filter metric data for this combination
          d_filt <- metric_dat %>%
            filter(HomeRangeSize == h, Threshold == th)
          
          if (nrow(d_filt) == 0) next
          
          # Use the yearly_aggregate_metric function from main script
          d_year <- tryCatch(
            yearly_aggregate_metric(d_filt, metric_col),
            error = function(e) {
              warning(sprintf("      Skipping %s: %s", metric_col, e$message))
              return(NULL)
            }
          )
          
          if (is.null(d_year) || nrow(d_year) == 0) next
          
          # Test each lag
          for (lag_y in lags) {
            
            # Lag the ice data
            d_lagged <- d_year %>%
              mutate(year = year + lag_y)
            
            # Merge with penguin data from BOTH regions
            combined_data <- penguin_full %>%
              left_join(d_lagged, by = "year") %>%
              filter(!is.na(metric_value), !is.na(growth_rate)) %>%
              mutate(region = factor(region))
            
            # Check if we have sufficient data
            if (nrow(combined_data) < min_rows_for_fit) next
            
            # Check that we have data from both regions
            regions_present <- combined_data %>%
              group_by(region) %>%
              summarise(n = n(), .groups = "drop")
            
            if (nrow(regions_present) < 2) next
            
            # Check for variation in metric
            if (sd(combined_data$metric_value, na.rm = TRUE) == 0) next
            
            # Fit interaction model
            tryCatch({
              
              # Full model with interaction
              model_interaction <- gls(
                growth_rate ~ metric_value * region,
                correlation = corAR1(form = ~1 | site_id),
                data = combined_data
              )
              
              # Extract coefficients
              coef_summary <- summary(model_interaction)$tTable
              
              # The interaction term tests if slopes differ between regions
              interaction_row <- grep("metric_value.*region|region.*metric_value", 
                                      rownames(coef_summary))
              
              if (length(interaction_row) > 0) {
                interaction_p <- coef_summary[interaction_row, "p-value"]
                interaction_coef <- coef_summary[interaction_row, "Value"]
                interaction_se <- coef_summary[interaction_row, "Std.Error"]
                interaction_t <- coef_summary[interaction_row, "t-value"]
              } else {
                interaction_p <- NA
                interaction_coef <- NA
                interaction_se <- NA
                interaction_t <- NA
              }
              
              # Main effect p-values
              metric_main_p <- coef_summary["metric_value", "p-value"]
              region_main_p <- if ("regionCentral_WAP" %in% rownames(coef_summary)) {
                coef_summary["regionCentral_WAP", "p-value"]
              } else {
                NA
              }
              
              # Store results
              result_df <- data.frame(
                File = basename(f),
                Metric = metric_name,
                MetricColumn = metric_col,
                HomeRangeSize = h,
                Threshold = th,
                Lag = lag_y,
                N_Obs_Total = nrow(combined_data),
                N_Obs_Bransfield = sum(combined_data$region == "Bransfield"),
                N_Obs_Central_WAP = sum(combined_data$region == "Central_WAP"),
                N_Colonies = n_distinct(combined_data$site_id),
                AIC = AIC(model_interaction),
                BIC = BIC(model_interaction),
                Metric_Main_Effect_p = metric_main_p,
                Region_Main_Effect_p = region_main_p,
                Interaction_Coefficient = interaction_coef,
                Interaction_StdError = interaction_se,
                Interaction_tValue = interaction_t,
                Interaction_pValue = interaction_p,
                stringsAsFactors = FALSE
              )
              
              interaction_results_list[[length(interaction_results_list) + 1]] <- result_df
              
            }, error = function(e) {
              cat(sprintf("    Error fitting model (HR=%s, Th=%s, Lag=%d): %s\n",
                          h, th, lag_y, e$message))
            })
          }
        }
      }
    }
  }
  
  # Compile results
  if (length(interaction_results_list) == 0) {
    cat("\nNo interaction models successfully fitted.\n")
    return(NULL)
  }
  
  all_interactions <- bind_rows(interaction_results_list)
  
  # Apply Bonferroni correction
  n_tests <- nrow(all_interactions)
  all_interactions <- all_interactions %>%
    mutate(
      Bonferroni_Adjusted_pValue = pmin(Interaction_pValue * n_tests, 1)
    )
  
  # Identify significant differences
  significant_interactions <- all_interactions %>%
    filter(Bonferroni_Adjusted_pValue < 0.05) %>%
    arrange(Bonferroni_Adjusted_pValue)
  
  # Save results
  write.csv(all_interactions,
            file.path(interaction_dir, "interaction_tests_all.csv"),
            row.names = FALSE)
  
  write.csv(significant_interactions,
            file.path(interaction_dir, "interaction_tests_significant.csv"),
            row.names = FALSE)
  
  # Print summary
  cat("\n", rep("=", 70), "\n")
  cat("INTERACTION ANALYSIS SUMMARY\n")
  cat(rep("=", 70), "\n\n")
  cat(sprintf("Total interaction tests performed: %d\n", n_tests))
  cat(sprintf("Bonferroni correction applied: n = %d\n", n_tests))
  cat(sprintf("Significant regional differences (p < 0.05): %d\n", 
              nrow(significant_interactions)))
  
  if (nrow(significant_interactions) > 0) {
    cat("\n--- Metrics with Significant Regional Differences ---\n\n")
    for (i in 1:min(10, nrow(significant_interactions))) {
      row <- significant_interactions[i, ]
      cat(sprintf("%d. %s (HR=%s, Th=%s, Lag=%d)\n",
                  i, row$Metric, row$HomeRangeSize, row$Threshold, row$Lag))
      cat(sprintf("   Interaction p = %.4f (Bonferroni adjusted p = %.4f)\n",
                  row$Interaction_pValue, row$Bonferroni_Adjusted_pValue))
      cat(sprintf("   Interaction coefficient = %.4f ± %.4f\n\n",
                  row$Interaction_Coefficient, row$Interaction_StdError))
    }
    
    if (nrow(significant_interactions) > 10) {
      cat(sprintf("   ... and %d more (see CSV output)\n", 
                  nrow(significant_interactions) - 10))
    }
  } else {
    cat("\nNo metrics showed significantly different effects between regions\n")
    cat("after Bonferroni correction.\n")
  }
  
  cat("\nResults saved to:\n")
  cat(sprintf("  - %s\n", file.path(interaction_dir, "interaction_tests_all.csv")))
  cat(sprintf("  - %s\n", file.path(interaction_dir, "interaction_tests_significant.csv")))
  cat("\n", rep("=", 70), "\n")
  
  return(all_interactions)
}