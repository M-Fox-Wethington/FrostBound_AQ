# Load necessary libraries
library(ggplot2)

# Plot histogram for AMSR dataset
ggplot(data, aes(x = AMSR)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of AMSR Data",
    x = "AMSR Values (0-1)",
    y = "Count"
  ) +
  theme_minimal()

# Plot histogram for NSIDC dataset
ggplot(data, aes(x = NSIDC)) +
  geom_histogram(binwidth = 0.05, fill = "darkorange", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of NSIDC Data",
    x = "NSIDC Values (0-1)",
    y = "Count"
  ) +
  theme_minimal()


# Find the unique values in the AMSR dataset
unique_amsr_values <- sort(unique(data$AMSR))

# Display the first few unique values starting from zero
head(unique_amsr_values)
