# Load necessary libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Create binned categories for AMSR and NSIDC
data$AMSR_bin <- cut(data$AMSR, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
data$NSIDC_bin <- cut(data$NSIDC, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)

# Calculate the absolute error for each data point
data$Abs_Error <- abs(data$AMSR - data$NSIDC)

# Aggregate data to calculate mean absolute error and total error for each bin combination
error_data <- data %>%
  group_by(AMSR_bin, NSIDC_bin) %>%
  summarize(
    Mean_Abs_Error = mean(Abs_Error),  # Mean absolute error in each bin
    Total_Error = sum(Abs_Error),      # Total error for all points in the bin
    Count = n()                        # Number of points in the bin
  ) %>%
  na.omit()

# Calculate the relative contribution of each bin to the total error
total_error_sum <- sum(error_data$Total_Error)
error_data <- error_data %>%
  mutate(Relative_Error_Contribution = (Total_Error / total_error_sum) * 100)

# Plot the heatmap of mean absolute error with Brewer Blues
ggplot(error_data, aes(x = AMSR_bin, y = NSIDC_bin, fill = Mean_Abs_Error)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", name = "Mean Abs Error") +
  labs(
    title = "Heatmap of Mean Absolute Error Between AMSR and NSIDC",
    x = "Binned AMSR (0-1)",
    y = "Binned NSIDC (0-1)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot the heatmap of relative error contribution with Brewer Blues
ggplot(error_data, aes(x = AMSR_bin, y = NSIDC_bin, fill = Relative_Error_Contribution)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", name = "Relative Error (%)") +
  labs(
    title = "Relative Error Contribution by Binned AMSR and NSIDC Values",
    x = "Binned AMSR (0-1)",
    y = "Binned NSIDC (0-1)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
