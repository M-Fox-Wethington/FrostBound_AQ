# Load necessary libraries
library(data.table)
library(dplyr)

# Define the directory paths
results_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/results/model-results"
significant_results_dir <- file.path(results_dir, "bonferroni_significant_results")

# Create the new subdirectory if it doesn't exist
dir.create(significant_results_dir, showWarnings = FALSE)

# Function to clean HomeRangeSize column
clean_home_range_size <- function(df) {
  df %>%
    mutate(HomeRangeSize = as.numeric(gsub("km", "", HomeRangeSize)))
}

# Function to create subsets of data based on Bonferroni significance criteria and Metric
subset_significant_results <- function(file_path) {
  # Load the data
  data <- fread(file_path)
  
  # Clean the HomeRangeSize column
  data <- clean_home_range_size(data)
  
  # Get unique metrics
  metrics <- unique(data$Metric)
  
  # Initialize a list to store subsets
  metric_data_list <- list()
  
  # Loop through each metric and create subset
  for (metric in metrics) {
    metric_data <- data %>%
      filter(Metric == metric, Bonferroni_Adjusted_pValue < 0.05)
    
    if (nrow(metric_data) > 0) {
      metric_data_list[[metric]] <- metric_data
    }
  }
  
  return(metric_data_list)
}

# Get a list of CSV files in the directory
csv_files <- list.files(path = results_dir, pattern = "*.csv", full.names = TRUE)

# Apply the function to each CSV file and bind the results by Metric
all_significant_results <- lapply(csv_files, subset_significant_results)
all_significant_results <- do.call(rbind, unlist(all_significant_results, recursive = FALSE))

# Split the data by Metric and save to separate CSV files
split_data <- split(all_significant_results, all_significant_results$Metric)

# Save each subset to a separate CSV file
lapply(names(split_data), function(metric) {
  output_file_name <- paste0(metric, "_bonferroni_significant.csv")
  output_file_path <- file.path(significant_results_dir, output_file_name)
  write.csv(split_data[[metric]], output_file_path, row.names = FALSE)
})
