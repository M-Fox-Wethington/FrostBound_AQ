# Load necessary libraries
library(dplyr)
library(readr)
library(ggplot2)

# Input: local copy of model results (avoids OneDrive sync issues)
input_dir <- "C:/Users/michael.wethington.BRILOON/Documents/GLS_model_results"

# Output directory
output_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Supplementary figure name mapping: key = "metric_threshold"
figure_map <- list(
  "mean_persistence_0.15" = list(num = "S06", label = "Persistence_15pct"),
  "mean_persistence_0.3"  = list(num = "S07", label = "Persistence_30pct"),
  "mean_persistence_0.5"  = list(num = "S08", label = "Persistence_50pct"),
  "mean_duration_0.15"    = list(num = "S09", label = "Duration_15pct"),
  "mean_duration_0.3"     = list(num = "S10", label = "Duration_30pct"),
  "mean_duration_0.5"     = list(num = "S11", label = "Duration_50pct"),
  "sd_persistence_0.15"   = list(num = "S14", label = "Persistence_Variability_15pct"),
  "sd_persistence_0.5"    = list(num = "S15", label = "Persistence_Variability_50pct"),
  "sd_duration_0.15"      = list(num = "S16", label = "Duration_Variability_15pct"),
  "sd_duration_0.3"       = list(num = "S17", label = "Duration_Variability_30pct"),
  "sd_duration_0.5"       = list(num = "S18", label = "Duration_Variability_50pct")
)

# Define the metrics to analyze
metrics <- c("mean_duration", "sd_duration", "mean_persistence", "sd_persistence")

# Load and merge CSVs
all_significant_results <- list()
for (metric in metrics) {
  file_path <- file.path(input_dir, paste0("significant_model_results_indiv_", metric, ".csv"))
  if (file.exists(file_path)) {
    data <- read_csv(file_path) %>% mutate(Metric_Type = metric)
    all_significant_results[[metric]] <- data
  } else {
    cat("File does not exist:", file_path, "\n")
  }
}
merged_significant_results <- bind_rows(all_significant_results)

# Save merged CSV
merged_file_path <- file.path(output_dir, "merged_significant_model_results_duration_persistence.csv")
write_csv(merged_significant_results, merged_file_path)
cat("Merged CSV saved to:", merged_file_path, "\n")

# Load merged data for plotting
data <- read_csv(merged_file_path)

# Recode factors
data$HomeRangeSize <- factor(data$HomeRangeSize, levels = c("25km", "50km", "100km", "150km", "200km", "250km", "300km", "350km", "400km", "450km", "500km"))
data$Metric <- recode(data$Metric, 
                      'mean_duration' = 'Duration', 
                      'sd_duration' = 'Duration SD', 
                      'mean_persistence' = 'Open Water Frequency',
                      'sd_persistence' = 'Open Water Frequency SD')
data$Lag <- recode(data$Lag, 
                   'Lag_Indiv 1' = '1 Year Lag', 
                   'Lag_Indiv 2' = '2 Year Lag', 
                   'Lag_Indiv 3' = '3 Year Lag', 
                   'Lag_Indiv 4' = '4 Year Lag',
                   'Lag_Indiv 5' = '5 Year Lag')

# Function to create and save forest plots
create_forest_plot <- function(data, fig_info, output_dir) {
  if (nrow(data) == 0) return(NULL)
  
  filename <- paste0("Figure_", fig_info$num, "_", fig_info$label)
  
  plot <- ggplot(data, aes(x = HomeRangeSize, y = Coefficient, color = Metric)) +
    geom_point() +
    geom_errorbar(aes(ymin = Coefficient - StdError, ymax = Coefficient + StdError), width = 0.2) +
    facet_wrap(~ Lag, scales = 'free_x', nrow = 2) +
    labs(x = "Home Range Size", y = "Coefficient") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  print(plot)
  
  ggsave(file.path(output_dir, paste0(filename, ".png")), plot, width = 10, height = 8)
  pdf(file.path(output_dir, paste0(filename, ".pdf")), width = 10, height = 8)
  print(plot)
  dev.off()
  
  cat("Saved:", filename, "\n")
}

# Generate forest plots using the figure map
thresholds <- unique(data$Threshold)
for (metric in metrics) {
  for (threshold in thresholds) {
    key <- paste0(metric, "_", threshold)
    fig_info <- figure_map[[key]]
    
    if (is.null(fig_info)) {
      cat("No figure mapping for:", key, "- skipping\n")
      next
    }
    
    metric_data <- data %>% filter(Metric_Type == metric, Threshold == threshold)
    
    if (nrow(metric_data) == 0) {
      cat("No significant data for:", key, "- skipping\n")
      next
    }
    
    create_forest_plot(metric_data, fig_info, output_dir)
  }
}

cat("\nFigures saved to:", output_dir, "\n")
