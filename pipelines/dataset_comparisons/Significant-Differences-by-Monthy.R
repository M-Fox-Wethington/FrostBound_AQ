library(terra)
library(lubridate)
library(RColorBrewer)  # for color palettes

# Assuming raster1 and raster2 are already loaded
# Set up the layout for a 2x2 grid of plots
par(mfrow=c(2, 2), mar=c(3, 3, 2, 1))

# Loop through each month from June to September for the year 2014
for (month in 6:9) {
  # Extract the subset of data for the given month and year
  time_data1 <- time(raster1)
  date_vector1 <- as.Date(time_data1)
  month_indices1 <- which(month(date_vector1) == month & year(date_vector1) == 2014)
  
  # If there is data for the month, proceed with calculation and plotting
  if (length(month_indices1) > 0) {
    month_mean1 <- mean(raster1[[month_indices1]], na.rm = TRUE)
    month_mean1_resampled <- resample(month_mean1, raster2, method="bilinear")
    
    # Repeat process for raster2
    time_data2 <- time(raster2)
    date_vector2 <- as.Date(time_data2)
    month_indices2 <- which(month(date_vector2) == month & year(date_vector2) == 2014)
    month_mean2 <- mean(raster2[[month_indices2]], na.rm = TRUE)
    
    # Calculate differences and find significant changes
    differences <- month_mean1_resampled - month_mean2
    significant_diff_values <- ifel(abs(differences) > 25, 1, NA)
    
    # Plot the significant changes
    plot_title <- sprintf("Significant Differences for %04d-%02d", 2014, month)
    plot(significant_diff_values, main=plot_title, col="red")
    
    # New Code: Plotting the over- and underestimation
    # Use a diverging color scale to indicate positive and negative differences
    colors <- brewer.pal(11, "RdBu")
    breaks <- seq(min(differences[], na.rm=TRUE), max(differences[], na.rm=TRUE), length.out = length(colors) + 1)
    colfunc <- colorRampPalette(colors)
    plot(differences, col=colfunc(length(colors)), breaks=breaks, main="Over/Underestimation", zlim=range(breaks))
  }
}

# Reset the plotting layout
par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
