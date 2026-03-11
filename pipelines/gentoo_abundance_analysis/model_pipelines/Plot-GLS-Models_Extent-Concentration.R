# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Input: local copy of model results (avoids OneDrive sync issues)
input_dir <- "C:/Users/michael.wethington.BRILOON/Documents/GLS_model_results/concentration_extent"

# Output directory
output_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Supplementary figure name mapping
figure_map <- list(
  "Sea Ice Extent"        = list(num = "S04", label = "Sea_Ice_Extent"),
  "Sea Ice Concentration" = list(num = "S05", label = "Sea_Ice_Concentration"),
  "Sea Ice Extent SD"     = list(num = "S12", label = "Sea_Ice_Extent_SD"),
  "Concentration SD"      = list(num = "S13", label = "Concentration_SD")
)

metrics <- names(figure_map)

# Load the data
data <- read_csv(file.path(input_dir, "significant_model_results_indiv_ice_metrics_extent_extentsd.csv"))

# Recode factors
data$HomeRangeSize <- factor(data$HomeRangeSize, levels = c("25km", "50km", "100km", "150km", "200km", "250km", "300km", "350km", "400km", "450km", "500km"))
data$Metric <- recode(data$Metric, 
                      'mean_sic' = 'Sea Ice Concentration', 
                      'sd_sic' = 'Concentration SD', 
                      'ice_extent_km2' = 'Sea Ice Extent',
                      'Extent_SD' = 'Sea Ice Extent SD')
data$Lag <- recode(data$Lag, 
                   'Lag_Indiv 1' = '1 Year Lag', 
                   'Lag_Indiv 2' = '2 Year Lag', 
                   'Lag_Indiv 3' = '3 Year Lag', 
                   'Lag_Indiv 4' = '4 Year Lag', 
                   'Lag_Indiv 5' = '5 Year Lag')

# Function to create and save forest plots
create_forest_plot <- function(data, metric, fig_info, output_dir) {
  if (nrow(data) == 0) {
    cat("No data available for", metric, "\n")
    return(NULL)
  }
  
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

# Generate forest plots
for (metric in metrics) {
  metric_data <- data %>% filter(Metric == metric)
  create_forest_plot(metric_data, metric, figure_map[[metric]], output_dir)
}

cat("\nFigures saved to:", output_dir, "\n")
