# Load necessary libraries
library(dplyr)
library(ggplot2)

# Load your Gentoo penguin dataset
gentoo_file_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/GentooCounts_48_1.csv"
gentoo_data <- read.csv(gentoo_file_path)

# Aggregate total Gentoo penguin abundance per year
annual_abundance <- gentoo_data %>%
  group_by(season_starting) %>%
  summarise(total_abundance = sum(penguin_count, na.rm = TRUE))

# Plot the total Gentoo penguin abundance over time
plot <- ggplot(annual_abundance, aes(x = season_starting, y = total_abundance)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Gentoo Penguin Abundance Over Time",
    x = "Year",
    y = "Total Gentoo Penguin Abundance"
  ) +
  theme_minimal()

# Save the plot
ggsave("D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo_abundance_trend.png", plot)

# Display the plot
print(plot)
