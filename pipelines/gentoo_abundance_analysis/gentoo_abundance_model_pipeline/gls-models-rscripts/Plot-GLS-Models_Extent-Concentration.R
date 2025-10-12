# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Define the directory containing the significant CSV files
results_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/results/model-results/concentration_extent"

# Define the metrics to analyze
metrics <- c("Sea Ice Concentration", "Concentration SD", "Sea Ice Extent", "Sea Ice Extent SD")

# Load the data
data <- read_csv(file.path(results_dir, "significant_model_results_indiv_ice_metrics_extent_extentsd.csv"))

# Reorder the factor levels for HomeRangeSize and Metric
data$HomeRangeSize <- factor(data$HomeRangeSize, levels = c("25km", "50km", "100km", "150km", "200km", "250km", "300km", "350km", "400km", "450km", "500km"))
data$Metric <- recode(data$Metric, 
                      'mean_sic' = 'Sea Ice Concentration', 
                      'sd_sic' = 'Concentration SD', 
                      'ice_extent_km2' = 'Sea Ice Extent',
                      'Extent_SD' = 'Sea Ice Extent SD')

# Recode the Lag column
data$Lag <- recode(data$Lag, 
                   'Lag_Indiv 1' = '1 Year Lag', 
                   'Lag_Indiv 2' = '2 Year Lag', 
                   'Lag_Indiv 3' = '3 Year Lag', 
                   'Lag_Indiv 4' = '4 Year Lag', 
                   'Lag_Indiv 5' = '5 Year Lag')

# Function to create and save forest plots
create_forest_plot <- function(data, metric, results_dir) {
  if (nrow(data) == 0) {
    cat("No data available for", metric, "\n")
    return(NULL)
  }
  
  plot <- ggplot(data, aes(x = HomeRangeSize, y = Coefficient, color = Metric)) +
    geom_point() +
    geom_errorbar(aes(ymin = Coefficient - StdError, ymax = Coefficient + StdError), width = 0.2) +
    facet_wrap(~ Lag, scales = 'free_x', nrow = 2) +
    labs(title = paste("Forest Plot of Coefficients for", metric),
         x = "Home Range Size",
         y = "Coefficient") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot
  print(plot)
  
  # Save the plot as PNG
  ggsave(file.path(results_dir, paste0("forest_plot_", metric, ".png")), plot, width = 10, height = 8)
  
  # Save the plot as PDF
  pdf(file.path(results_dir, paste0("forest_plot_", metric, ".pdf")), width = 10, height = 8)
  print(plot)
  dev.off()
}

# Generate forest plots for each metric
for (metric in metrics) {
  print(metric)
  metric_data <- data %>% filter(Metric == metric)
  create_forest_plot(metric_data, metric, results_dir)
}

# Generate the combined forest plot
combined_plot <- ggplot(data, aes(x = HomeRangeSize, y = Coefficient, color = Metric)) +
  geom_point() +
  geom_errorbar(aes(ymin = Coefficient - StdError, ymax = Coefficient + StdError), width = 0.2) +
  facet_wrap(~ Lag, scales = 'free_x', nrow = 2) +
  labs(title = "Forest Plot of Coefficients by Metric, Home Range Size, and Lag",
       x = "Home Range Size",
       y = "Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the combined plot
print(combined_plot)

# Save the combined plot as PNG
ggsave(file.path(results_dir, "forest_plot_combined_metrics_Concentration-Extent.png"), combined_plot, width = 12, height = 10)

# Save the combined plot as PDF
pdf(file.path(results_dir, "forest_plot_combined_metrics_Concentration-Extent.pdf"), width = 12, height = 10)
print(combined_plot)
dev.off()
