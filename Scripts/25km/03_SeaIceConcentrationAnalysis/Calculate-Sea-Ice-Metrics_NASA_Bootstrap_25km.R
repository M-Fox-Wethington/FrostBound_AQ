library(terra)
library(doParallel)

# Define the file path for the raster stack
r_fp <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/NASA_Bootstrap_sic/raster-stack/NASA_Bootstrap_SIC_Time_Series.nc"

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
    
    # Calculate sea ice extent
    # Convert sea ice concentration to binary (1 if >= 0.15, 0 otherwise)
    binary_ice_values <- ifelse(data_values >= 0.15, 1, 0)
    binary_ice_raster <- layer_data
    values(binary_ice_raster) <- binary_ice_values
    
    # Calculate total ice extent by summing the areas of cells marked as ice
    cell_area = 156.25  # Area of each cell in square kilometers
    ice_extent <- sum(values(binary_ice_raster) * cell_area, na.rm = TRUE)
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
