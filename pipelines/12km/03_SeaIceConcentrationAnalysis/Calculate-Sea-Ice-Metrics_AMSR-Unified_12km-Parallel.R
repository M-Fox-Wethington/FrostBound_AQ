library(terra)        # For raster data manipulation
library(doParallel)   # For parallel processing
library(foreach)      # For looping with parallel support

# Function to calculate mean SIC and sea ice extent for a given raster layer within a file
calculate_layer_stats <- function(layer_index, file_path) {
  raster_stack <- rast(file_path)                  # Dynamically load the raster stack from file
  layer_data <- raster_stack[[layer_index]]        # Extract the specific layer
  layer_date <- time(layer_data)                   # Extract date from layer metadata
  
  # Replace 0 with NA
  values(layer_data)[values(layer_data) == 0] <- NA
  
  # Calculate mean SIC excluding NAs and determine ice extent
  mean_sic <- as.numeric(terra::global(layer_data, fun = 'mean', na.rm = TRUE))
  cell_area_sq_meters <- prod(res(layer_data))    # Calculate the area of one cell in square meters
  valid_ice_cells <- sum(values(layer_data) >= 15, na.rm = TRUE)  # Count cells with valid ice concentration
  
  # Convert total ice area to square kilometers
  total_ice_area_sq_km <- (valid_ice_cells * cell_area_sq_meters) / 1e6
  
  rm(raster_stack, layer_data)   # Clear memory
  gc()                           # Run garbage collection
  
  return(list(mean_sic = mean_sic, ice_extent_km = total_ice_area_sq_km, date = layer_date))
}

# Directory and file path setup
r_fp <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/AMSR-Unified_SIC_Full-Catalog_Time_Series.nc"

# Set up parallel computing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Send necessary function and file path to all workers
clusterExport(cl, varlist = c("calculate_layer_stats", "r_fp"))

# Perform calculations in parallel with error handling
results <- foreach(i = 1:nlyr(rast(r_fp)), .packages = 'terra', .errorhandling = "pass") %dopar% {
  tryCatch({
    calculate_layer_stats(i, r_fp)
  }, error = function(e) {
    return(list(error = TRUE, message = e$message))
  })
}

# Stop the cluster after the task is done to free up resources
stopCluster(cl)

# Handling results and constructing the final dataframe
result_df <- do.call(rbind, lapply(results, function(x) {
  if (!is.null(x$error)) {
    return(data.frame(Date = NA, MeanSIC = NA, IceExtent_km = NA, error = x$message))
  } else {
    return(data.frame(Date = x$date, MeanSIC = x$mean_sic, IceExtent_km = x$ice_extent_km, error = NA))
  }
}))

# Save results to CSV and RDS for further use
write.csv(result_df, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/12km_AMSR-Unified_Metrics.csv", row.names = FALSE)
saveRDS(result_df, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/12km_AMSR-Unified_Metrics.rds")

# Display some of the results for verification
print(head(result_df))
