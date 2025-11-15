# Load necessary libraries
library(terra)
library(sf)
library(dplyr)
library(nlme)
library(lubridate)
library(data.table)

# General function to compute annual ice metrics
compute_annual_ice_metrics <- function(ice_data, metrics) {
  ice_data <- ice_data %>%
    mutate(date = as.Date(date, format = "%m/%d/%Y"),
           year = year(date))
  
  annual_ice_metrics <- ice_data %>%
    group_by(year, Threshold, HomeRangeSize) %>%
    summarize(across(all_of(metrics), ~mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
              across(all_of(metrics), ~sd(.x, na.rm = TRUE), .names = "sd_{.col}"),
              .groups = 'drop')
  
  return(annual_ice_metrics)
}

# Function to get individual precomputed ice metrics
get_individual_ice_metrics <- function(year, ice_data, metric) {
  ice_value <- ice_data %>%
    filter(year == !!year) %>%
    pull(!!sym(metric))
  
  if (length(ice_value) > 0) {
    return(ice_value)
  } else {
    return(NA)
  }
}

# Function to fit GLS models for individual lag years
fit_gls_models_indiv <- function(penguin_data, ice_data, metric, threshold, home_range_size) {
  results <- list()
  
  for (lag in 1:5) {
    penguin_data <- penguin_data %>%
      rowwise() %>%
      mutate(overwinter_ice = get_individual_ice_metrics(year - lag, ice_data, metric)) %>%
      ungroup() %>%
      filter(!is.na(overwinter_ice) & !is.na(growth_rate) & growth_rate <= 3)
    
    if (nrow(penguin_data) < 3) next
    
    print(paste("Fitting GLS model for individual lag", lag, "year(s)"))
    formula <- as.formula("growth_rate ~ overwinter_ice")
    model <- gls(formula, correlation = corAR1(form = ~1 | site_id), data = penguin_data)
    summary_model <- summary(model)
    bonferroni_p_value <- p.adjust(summary_model$tTable[2, 4], method = "bonferroni", n = 5)
    
    if (bonferroni_p_value < 0.05) {
      results[[paste("Lag_Indiv", lag)]] <- list(
        AIC = AIC(model),
        BIC = BIC(model),
        coefficients = summary_model$tTable,
        p_value = summary_model$tTable[2, 4],
        bonferroni_p_value = bonferroni_p_value,
        model = model,
        penguin_data = penguin_data
      )
    }
  }
  
  return(results)
}

# Main function to run GLS analysis for each condition and save results
run_gls_analysis <- function(penguin_abundance_data, ice_data, metrics, results_dir) {
  # Initialize results list
  all_gls_results_indiv <- list()
  all_results_df_indiv <- list()
  significant_results_df_indiv <- list()
  
  # Compute annual ice metrics
  annual_ice_metrics <- compute_annual_ice_metrics(ice_data, metrics)
  
  # List unique thresholds and home range sizes
  thresholds <- unique(annual_ice_metrics$Threshold)
  home_range_sizes <- unique(annual_ice_metrics$HomeRangeSize)
  
  for (metric in metrics) {
    for (threshold in thresholds) {
      for (home_range_size in home_range_sizes) {
        ice_subset <- annual_ice_metrics %>%
          filter(Threshold == threshold, HomeRangeSize == home_range_size)
        
        print(paste("Running GLS analysis for metric:", metric, "threshold:", threshold, "home range size:", home_range_size))
        
        # Perform GLS analysis for individual lag years
        gls_results_indiv <- fit_gls_models_indiv(penguin_abundance_data, ice_subset, paste("mean_", metric, sep = ""), threshold, home_range_size)
        
        # Store results for individual lag years
        all_gls_results_indiv[[paste("Annual_Indiv", metric, threshold, home_range_size, sep = "_")]] <- gls_results_indiv
        
        results_list_indiv <- lapply(names(gls_results_indiv), function(lag) {
          result <- gls_results_indiv[[lag]]
          data.frame(
            Metric = metric,
            HomeRangeSize = home_range_size,
            Threshold = threshold,
            Lag = lag,
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
            Bonferroni_Adjusted_pValue = result$bonferroni_p_value
          )
        })
        
        all_results_df_indiv <- append(all_results_df_indiv, results_list_indiv)
        significant_results_df_indiv <- append(significant_results_df_indiv, Filter(function(x) x$Bonferroni_Adjusted_pValue < 0.05, results_list_indiv))
      }
    }
  }
  
  # Export the compiled results to CSV
  all_results_df_indiv <- do.call(rbind, all_results_df_indiv)
  significant_results_df_indiv <- do.call(rbind, significant_results_df_indiv)
  
  write.csv(all_results_df_indiv, file.path(results_dir, "model_results_indiv_ice_metrics.csv"), row.names = FALSE)
  write.csv(significant_results_df_indiv, file.path(results_dir, "significant_model_results_indiv_ice_metrics.csv"), row.names = FALSE)
  
  # Save only significant results to RDS
  significant_gls_results_indiv <- list()
  for (metric in names(all_gls_results_indiv)) {
    significant_gls_results_indiv[[metric]] <- Filter(Negate(is.null), all_gls_results_indiv[[metric]])
  }
  saveRDS(significant_gls_results_indiv, file.path(results_dir, "gls_analysis_results_indiv.rds"))
}

### Load Data and Run Analysis

# Load Gentoo penguin data
penguin_abundance_data <- fread("D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/inputs/modeled_gentoo_parameters.csv")

str(penguin_abundance_data)

# Adjust the season column in the penguin abundance data
penguin_abundance_data <- penguin_abundance_data %>%
  mutate(year = 1970 + season - 1)

# Load the SIC statistics
metric_calculation_csv <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/metric-calculation-csv"
ice_data <- fread(file.path(metric_calculation_csv, "daily_sic_statistics.csv"))
str(ice_data)


# Define the metrics to analyze
metrics <- c("mean_sic", "sd_sic", "ice_extent_km2")

# Define the results directory
results_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/results/model-results"
dir.create(results_dir, showWarnings = FALSE)

# Run the GLS analysis on the full dataset
results_output <- run_gls_analysis(penguin_abundance_data, ice_data, metrics, results_dir)
