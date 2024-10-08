# Load necessary libraries
library(ggplot2)
library(scales)  # For rescale function

# Load the data
data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/trend-analysis/Trends/significant_results_with_CIs.csv")

# Ensure that relevant columns are not missing
data <- na.omit(data)

# Ensure the Metric and Region columns are correctly ordered and labeled
data$Metric <- factor(data$Metric, levels = c("Extent", "Extent (Var.)", 
                                              "SIC", "SIC (Var.)", 
                                              "Duration", "Duration (Var.)", 
                                              "Open Water Frequency", "Open Water Frequency (Var.)"))

data$Region <- factor(data$Region, levels = c("Complete Study Area", "Offshore", "Northern Shelf", "Middle Shelf", "Southern Shelf"))

# Replace "Annual" with "Winter" in the Month column
data$Month <- as.character(data$Month)
data$Month[data$Month == "Annual"] <- "Winter"

# Convert the Month column to a factor with custom labels
month_labels <- c("6" = "June", "7" = "July", "8" = "August", "9" = "September", "Winter" = "Winter")
data$Month <- factor(data$Month, levels = names(month_labels), labels = month_labels)

# Function to scale values while preserving the sign
relative_scale <- function(x) {
  min_x <- min(x)
  max_x <- max(x)
  
  # If all values are the same, return them as is
  if (min_x == max_x) {
    return(rep(0, length(x)))
  } else {
    scaled_x <- rescale(abs(x), to = c(0, 1))  # Scale absolute values between 0 and 1
    return(scaled_x * sign(x))  # Reapply the sign to the scaled values
  }
}

# Apply the relative scaling function to each Metric
data$Relative_Scaled_Slope <- ave(data$Time_Estimate, data$Metric, FUN = relative_scale)

# Check the summary of the relative scaled values
cat("\nSummary of Relative Scaled Slope:\n")
print(summary(data$Relative_Scaled_Slope))

# Inspect SIC-specific data to ensure it's correctly scaled
cat("\nSIC-specific Data Inspection:\n")
sic_data <- subset(data, Metric == "SIC")
print(sic_data[, c("Time_Estimate", "Relative_Scaled_Slope")])

# Create the heatmap with the relative scaled slopes and the custom color palette
ggplot(data, aes(x = Month, y = Region, fill = Relative_Scaled_Slope)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("#E46726", "white", "#6D9EC1"), 
                       name = "Trend Strength & Direction",
                       limits = c(-1, 1),  # Ensure color mapping from -1 to 1
                       oob = scales::squish) +
  labs(title = "Seasonal Trend Directions and Strength by Region", fill = "Relative Scaled Slope") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(2, "lines")) +
  facet_wrap(~ Metric, ncol = 2) +
  coord_fixed(ratio = 1)  # Set aspect ratio to make tiles more square
