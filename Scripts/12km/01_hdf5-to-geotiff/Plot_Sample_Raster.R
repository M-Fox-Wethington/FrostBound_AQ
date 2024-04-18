# Specify the path to the raster file
raster_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-Data/Processed/NSIDC-AU_SI12_SeaIce_12km_20120703.tif"  # Change YYYYMMDD to the actual date in your filename

# Load the raster file
sic_raster <- rast(raster_path)


# Check the raster and print basic information
print(sic_raster)

# Correct way to find the minimum and maximum values in the raster
raster_min <- min(values(sic_raster), na.rm = TRUE)
raster_max <- max(values(sic_raster), na.rm = TRUE)

# Output the min and max values using cat
cat("Minimum value in raster: ", raster_min, "\n")
cat("Maximum value in raster: ", raster_max, "\n")


# Plot the raster
# Ensure that raster_min and raster_max are not NA and are finite
if (!is.na(raster_min) && !is.na(raster_max) && is.finite(raster_min) && is.finite(raster_max)) {
  plot(sic_raster, 
       main = "Sea Ice Concentration", 
       col = hcl.colors(100, "Blues"),  # Using a color scale that is visually intuitive for water/ice
       breaks = seq(raster_min, raster_max, length.out = 101))  # Define breaks in the data for coloring
} else {
  cat("Unable to plot: non-finite min or max values.\n")
}