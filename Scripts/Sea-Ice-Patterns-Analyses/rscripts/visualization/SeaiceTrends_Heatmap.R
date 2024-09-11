# Load necessary libraries
library(ggplot2)
library(RColorBrewer)

# Load the data
data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/RStudioProject/Scripts/Sea-Ice-Patterns-Analyses/sea_ice_trends_table4.csv")

# View the first few rows to confirm it loaded correctly
head(data)

# Check the column names
colnames(data)

# Z-Score Standardization of Slopes within each Metric to preserve direction and magnitude
data$Standardized_Slope <- ave(data$`Slope..Time.`, data$Metric, 
                               FUN = function(x) scale(x))

# Reorder the regions on the y-axis
data$Region <- factor(data$Region, levels = c("Complete Area", "Offshore", "Northern Shelf", "Middle Shelf", "Southern Shelf"))

# Ensure related metrics are ordered correctly
data$Metric <- factor(data$Metric, levels = c("Extent", "Extent Variability", 
                                              "SIC", "SIC Variability", 
                                              "Duration", "Duration Variability", 
                                              "Open Water Frequency", "Open Water Frequency Var."))

# Create the heatmap with faceting, standardized slopes, increased spacing, and custom RdBu gradient
ggplot(data, aes(x = Month, y = Region, fill = Standardized_Slope)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = colorRampPalette(brewer.pal(3, "RdBu"))(256), 
                       name = "Standardized Slope (Trend Strength & Direction)") +
  labs(title = "Seasonal Trend Directions and Strength by Region ", fill = "Standardized Slope") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(2, "lines")) +  # Increase spacing between columns
  facet_wrap(~ Metric, ncol = 2)
