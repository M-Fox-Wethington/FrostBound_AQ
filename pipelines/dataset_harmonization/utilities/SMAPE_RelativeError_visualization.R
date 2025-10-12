# Install and load the necessary package for GSEA palette
# install.packages("ggsci")  # Uncomment this line if you don't have ggsci installed
library(ggsci)
library(ggplot2)
library(dplyr)


# Load the harmonized NSIDC dataset
cat("Loading the harmonized NSIDC dataset...\n")
harmonized_stack <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_harmonized_common_dates.tif")

# Load the AMSR dataset
cat("Loading the AMSR dataset...\n")
amsr_overlap <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/amsr_common_dates.tif")

# Ensure that the time information is consistent
time(harmonized_stack) <- time(amsr_overlap)

# Recast any values above 1 to 1
harmonized_stack <- clamp(harmonized_stack, lower = 0, upper = 1, values = TRUE)

# Prepare the data for plotting
data <- data.frame(
  AMSR = as.vector(values(amsr_overlap)),
  NSIDC = as.vector(values(harmonized_stack))
)
data <- na.omit(data)  # Remove rows with NA values




# Create binned categories for AMSR and NSIDC
data$AMSR_bin <- cut(data$AMSR, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
data$NSIDC_bin <- cut(data$NSIDC, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)

# Function to calculate SMAPE
calc_smape <- function(actual, predicted) {
  return(100 * mean(abs(actual - predicted) / ((abs(actual) + abs(predicted)) / 2)))
}

# Aggregate data to calculate mean SMAPE for each bin combination
error_data <- data %>%
  group_by(AMSR_bin, NSIDC_bin) %>%
  summarize(
    Mean_SMAPE = calc_smape(AMSR, NSIDC),  # Mean SMAPE for each bin
    Count = n()                            # Number of points in the bin
  ) %>%
  na.omit()

# Calculate the relative contribution of each bin to the total SMAPE
total_smape_sum <- sum(error_data$Mean_SMAPE * error_data$Count)
error_data <- error_data %>%
  mutate(Relative_Error_Contribution = (Mean_SMAPE * Count / total_smape_sum) * 100)

# Create the plot with the GSEA palette and error labels
ggplot(error_data, aes(x = AMSR_bin, y = NSIDC_bin, fill = Relative_Error_Contribution)) +
  geom_tile() +
  scale_fill_gsea(name = "Relative Error (%)") +
  labs(
    title = "Relative Error Contribution by Binned AMSR and NSIDC Values",
    x = "Binned AMSR (0-1)",
    y = "Binned NSIDC (0-1)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(label = round(Relative_Error_Contribution, 1)), color = "black", size = 3)  # Add error labels in each cell
