library(terra)
library(doParallel)

# Define the file path for the raster stack
r_fp <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/AMSR-Unified_SIC_Full-Catalog_Time_Series.nc"

# Load the raster stack outside of the parallel loop for main session setup
raster_stack <- rast(r_fp)

# Modify the function to work independently in each worker
calculate_layer_stats <- function(index, file_path) {
  local_raster <- rast(file_path)
  layer_data <- local_raster[[index]]
  
  # Cast values of 0 to NA
  values(layer_data)[values(layer_data) == 0] <- NA
  
  # Check if all values are NA
  mean_value <- NA
  ice_extent <- NA
  if (any(!is.na(values(layer_data)))) {
    mean_value <- mean(values(layer_data), na.rm = TRUE)
    binary_ice_values <- ifelse(values(layer_data) >= 0.15, 1, 0)
    binary_ice_raster <- layer_data
    values(binary_ice_raster) <- binary_ice_values
    ice_extent <- sum(values(binary_ice_raster) * res(binary_ice_raster)[1] * res(binary_ice_raster)[2], na.rm = TRUE)
  }
  
  return(list(mean = mean_value, ice_extent = ice_extent))
}

# Set up parallel computing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Send necessary function and file path to all workers
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
write.csv(result_df, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/12km_SIC.csv", row.names = FALSE)

# Save results to RDS
saveRDS(result_df, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/12km_SIC.rds")

# Display some of the results for verification
print(head(result_df))
