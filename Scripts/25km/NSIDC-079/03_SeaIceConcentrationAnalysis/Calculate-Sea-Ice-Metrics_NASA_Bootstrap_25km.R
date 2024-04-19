library(terra)
library(doParallel)

# Define the file path for the raster stack
r_fp <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/25km_NASA_Bootstrap/stack/NSIDC079_25km_SIC_Time_Series.nc"

# Load the raster stack outside of the parallel loop for main session setup
raster_stack <- rast(r_fp)

# Define your function to calculate mean and check for sea ice extent
calculate_layer_stats <- function(layer_index, file_path) {
  local_raster <- rast(file_path)  # Load raster within the function
  layer_data <- local_raster[[layer_index]]
  
  # Check if all values are NA
  data_values <- values(layer_data)
  mean_value <- NA
  ice_extent <- NA
  
  if (!all(is.na(data_values))) {
    # Calculate mean excluding NA values
    mean_value <- mean(data_values, na.rm = TRUE)
    values(sic_raster)[values(sic_raster) == 1200] <- NA  # Land or ice-shelf
    
    # Calculate sea ice extent
    # Convert sea ice concentration to binary (1 if >= 0.15, 0 otherwise)
    binary_ice_values <- ifelse(data_values >= 0.15, 1, 0)
    binary_ice_raster <- layer_data
    values(binary_ice_raster) <- binary_ice_values
    
    
    
    # Calculate total ice extent by summing the areas of cells marked as ice
    ice_extent <- terra::expanse(binary_ice_raster, unit = "m")
  }
  
  return(list(mean = mean_value, ice_extent = ice_extent))
}

# Set up parallel computing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
clusterExport(cl, varlist = c("calculate_layer_stats", "r_fp"))

# Extract dates for each layer
dates <- as.Date(time(raster_stack))

# Perform calculations in parallel
results <- foreach(i = 1:nlyr(raster_stack), .packages = 'terra') %dopar% {
  calculate_layer_stats(i, r_fp)
}

# Stop the cluster after the task is done to free up resources
stopCluster(cl)

# Create a dataframe with results and additional data
result_df <- data.frame(
  Index = 1:nlyr(raster_stack),
  Date = dates,
  MeanSIC = sapply(results, `[[`, "mean"),
  IceExtent = sapply(results, `[[`, "ice_extent")
)

# Save results to CSV
write.csv(result_df, "sea_ice_data.csv", row.names = FALSE)

# Save results to RDS
saveRDS(result_df, "sea_ice_data.rds")

# Display some of the results for verification
print(head(result_df))

# Assuming result_df is already loaded or computed

# Convert IceExtent from square meters to square kilometers
result_df$IceExtent_km2 <- result_df$IceExtent / 1e6

# Print the updated dataframe with IceExtent in square kilometers
print(head(result_df))


# Assuming that you have finished all your calculations and possibly closed your parallel computing setup
# Now, let's plot the first layer of the raster stack
layer_to_plot <- 1  # Change this to whatever layer you want to inspect

# Extract the specific layer from the raster stack
single_layer <- raster_stack[[layer_to_plot]]

# Plot the layer using terra's plot function
plot(single_layer, main = paste("Sea Ice Concentration Layer", layer_to_plot, format(as.Date(time(single_layer)), "%Y-%m-%d")))

# You might also want to add additional aesthetic settings or use a color palette that is more suitable for visualizing sea ice concentration
library(RColorBrewer)
col_breaks <- seq(0, 1, length.out = 11)  # Change as needed based on the data range
colors <- colorRampPalette(brewer.pal(9, "Blues"))(10)  # Choosing a color scheme

plot(single_layer, breaks = col_breaks, col = colors, main = paste("Sea Ice Concentration Layer", layer_to_plot, format(as.Date(time(single_layer)), "%Y-%m-%d")))

