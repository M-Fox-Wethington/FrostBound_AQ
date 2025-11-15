# Load necessary libraries
library(dplyr)
library(readr)
library(ggplot2)

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

# Load the merged data for plotting
data <- read_csv(merged_file_path)

# Reorder the factor levels for HomeRangeSize and Metric
data$HomeRangeSize <- factor(data$HomeRangeSize, levels = c("25km", "50km", "100km", "150km", "200km", "250km", "300km", "350km", "400km", "450km", "500km"))
data$Metric <- recode(data$Metric, 
                      'mean_duration' = 'Duration', 
                      'sd_duration' = 'Duration SD', 
                      'mean_persistence' = 'Open Water Frequency',
                      'sd_persistence' = 'Open Water Frequency SD')

# Recode the Lag column
data$Lag <- recode(data$Lag, 
                   'Lag_Indiv 1' = '1 Year Lag', 
                   'Lag_Indiv 2' = '2 Year Lag', 
                   'Lag_Indiv 3' = '3 Year Lag', 
                   'Lag_Indiv 4' = '4 Year Lag',
                   'Lag_Indiv 5' = '5 Year Lag')

# Function to create and save forest plots
create_forest_plot <- function(data, metric, threshold, results_dir) {
  if (nrow(data) == 0) {
    cat("No data available for", metric, "at", threshold * 100, "% threshold\n")
    return(NULL)
  }
  
  plot <- ggplot(data, aes(x = HomeRangeSize, y = Coefficient, color = Metric)) +
    geom_point() +
    geom_errorbar(aes(ymin = Coefficient - StdError, ymax = Coefficient + StdError), width = 0.2) +
    facet_wrap(~ Lag, scales = 'free_x', nrow = 2) +
    labs(title = paste("Forest Plot of Coefficients for", metric, "at", threshold * 100, "% Threshold"),
         x = "Home Range Size",
         y = "Coefficient") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot
  print(plot)
  
  # Save the plot as PNG
  ggsave(file.path(results_dir, paste0("forest_plot_", metric, "_", threshold * 100, "pct.png")), plot, width = 10, height = 8)
  
  # Save the plot as PDF
  pdf(file.path(results_dir, paste0("forest_plot_", metric, "_", threshold * 100, "pct.pdf")), width = 10, height = 8)
  print(plot)
  dev.off()
}

# Generate forest plots for each metric and threshold
thresholds <- unique(data$Threshold)
for (metric in metrics) {
  for (threshold in thresholds) {
    metric_data <- data %>% filter(Metric_Type == metric, Threshold == threshold)
    create_forest_plot(metric_data, metric, threshold, results_dir)
  }
}

# Generate the combined forest plot
combined_plot <- ggplot(data, aes(x = HomeRangeSize, y = Coefficient, color = Metric)) +
  geom_point() +
  geom_errorbar(aes(ymin = Coefficient - StdError, ymax = Coefficient + StdError), width = 0.2) +
  facet_wrap(~ Lag + Threshold, scales = 'free_x', nrow = 2) +
  labs(title = "Forest Plot of Coefficients by Metric, Home Range Size, Lag, and Threshold",
       x = "Home Range Size",
       y = "Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the combined plot
print(combined_plot)

# Save the combined plot as PNG
ggsave(file.path(results_dir, "Figure_forest_plot_combined_metrics_thresholds_Durations-Persistence.png"), combined_plot, width = 12, height = 10)

# Save the combined plot as PDF
pdf(file.path(results_dir, "Figure_forest_plot_combined_metrics_thresholds_Durations-Persistence.pdf"), width = 12, height = 10)
print(combined_plot)
dev.off()
