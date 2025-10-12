# Load necessary libraries
library(terra)
library(sf)
library(dplyr)
library(nlme)
library(lubridate)
library(data.table)

# Function to get individual precomputed ice metrics
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

# Function: Fit GLS models for individual lag years
fit_gls_models_indiv <- function(penguin_data, ice_data, metric) {
  results <- list()
  
  for (lag in 1:5) {
    penguin_data <- penguin_data %>%
      rowwise() %>%
      mutate(overwinter_ice = get_individual_ice_metrics(year - lag, ice_data, metric)) %>%
      ungroup() %>%
      filter(!is.na(overwinter_ice))
    
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
  
  # List unique thresholds and home range sizes
  thresholds <- unique(ice_data$Threshold)
  home_range_sizes <- unique(ice_data$HomeRangeSize)
  
  for (metric in metrics) {
    for (threshold in thresholds) {
      for (home_range_size in home_range_sizes) {
        ice_subset <- ice_data %>%
          filter(Threshold == threshold, HomeRangeSize == home_range_size)
        
        print(paste("Running GLS analysis for metric:", metric, "threshold:", threshold, "home range size:", home_range_size))
        
        # Perform GLS analysis for individual lag years
        gls_results_indiv <- fit_gls_models_indiv(penguin_abundance_data, ice_subset, metric)
        
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
  for (metric in metrics) {
    metric_results_df_indiv <- do.call(rbind, lapply(all_results_df_indiv, function(x) if (x$Metric == metric) x else NULL))
    metric_significant_results_df_indiv <- do.call(rbind, lapply(significant_results_df_indiv, function(x) if (x$Metric == metric) x else NULL))
    
    write.csv(metric_results_df_indiv, file.path(results_dir, paste0("model_results_indiv_", metric, ".csv")), row.names = FALSE)
    write.csv(metric_significant_results_df_indiv, file.path(results_dir, paste0("significant_model_results_indiv_", metric, ".csv")), row.names = FALSE)
  }
  
  # Save only significant results to RDS
  significant_gls_results_indiv <- list()
  for (metric in names(all_gls_results_indiv)) {
    significant_gls_results_indiv[[metric]] <- Filter(Negate(is.null), all_gls_results_indiv[[metric]])
  }
  saveRDS(significant_gls_results_indiv, file.path(results_dir, "gls_analysis_results_indiv.rds"))
}

# Example Usage for a Small Subset

# Load Gentoo penguin data
penguin_abundance_data <- fread("D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/inputs/modeled_gentoo_parameters.csv")

# Adjust the season column in the penguin abundance data
penguin_abundance_data <- penguin_abundance_data %>%
  mutate(year = 1970 + season - 1)

# Load the SIC statistics
metric_calculation_csv <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/metric-calculation-csv"
ice_data <- fread(file.path(metric_calculation_csv, "sea_ice_duration_persistence_stats.csv"))

# Define the metrics to analyze
metrics <- c("mean_duration", "sd_duration", "mean_persistence", "sd_persistence")

# Define the results directory
results_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/results/model-results"
dir.create(results_dir, showWarnings = FALSE)

run_gls_analysis(penguin_abundance_data, ice_data, metrics, results_dir)




# Load necessary libraries
library(dplyr)
library(readr)

# Define the directory containing the significant CSV files
results_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/results/model-results"

# Define the metrics to analyze
metrics <- c("mean_duration", "sd_duration", "mean_persistence", "sd_persistence")

# Initialize an empty list to store the data frames
all_significant_results <- list()

# Loop through each metric and read the corresponding CSV file
for (metric in metrics) {
  file_path <- file.path(results_dir, paste0("significant_model_results_indiv_", metric, ".csv"))
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the CSV file and add the metric as a new column
    data <- read_csv(file_path) %>%
      mutate(Metric_Type = metric)
    all_significant_results[[metric]] <- data
  } else {
    cat("File does not exist:", file_path, "\n")
  }
}

# Combine all data frames into a single data frame
merged_significant_results <- bind_rows(all_significant_results)

# Save the merged data frame to a new CSV file
merged_file_path <- file.path(results_dir, "merged_significant_model_results_duration_persistence.csv")
write_csv(merged_significant_results, merged_file_path)

cat("Merged CSV file saved to:", merged_file_path, "\n")

