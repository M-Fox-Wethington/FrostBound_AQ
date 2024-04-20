library(terra)
library(doParallel)
library(foreach)

# Function to calculate mean SIC and sea ice extent
# This function includes detailed documentation for each step.
calculate_layer_stats <- function(layer_index, file_path) {
  # Load the specific layer from the raster stack to ensure efficient memory use
  raster_stack <- rast(file_path)
  layer_data <- raster_stack[[layer_index]]
  
  # Extract the date associated with the layer
  layer_date <- time(layer_data)
  
  # Calculate mean SIC, considering only values that are not NA
  mean_sic <- as.numeric(terra::global(layer_data, fun = 'mean', na.rm = TRUE))
  
  # Calculate the total area covered by valid sea ice cells
  # First, calculate the cell resolution (in meters) of input raster product, which is the area of one cell
  cell_area_sq_meters <- prod(res(layer_data))
  
  # Count the number of valid ice concentration cells (above 15%)
  valid_ice_cells <- sum(values(layer_data) >= 15.0, na.rm = TRUE)
  
  # Calculate the total area covered by valid sea ice cells in square kilometers
  total_ice_area_sq_km <- (valid_ice_cells * cell_area_sq_meters) / 1e6
  
  # Explicitly clear memory of loaded raster data
  rm(raster_stack, layer_data)
  gc()  # Force garbage collection
  
  return(list(mean_sic = mean_sic, ice_extent_km = total_ice_area_sq_km, date = layer_date))
}

# File path for the raster stack
r_fp <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/tm/Chunk_1978-10-26.nc"

# Set up parallel computing environment
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Export necessary variables and functions to each worker
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

# Compile results into a data frame if no errors
result_df <- data.frame(
  Index = 1:length(results),
  Date = sapply(results, function(x) as.character(x$date)),  # Convert date to character for dataframe
  MeanSIC = sapply(results, function(x) x$mean_sic),
  IceExtent_km = sapply(results, function(x) x$ice_extent_km)
)

# Save results to CSV
write.csv(result_df, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/tm/sea_ice_data.csv", row.names = FALSE)

# Optionally, load and verify one of the outputs
if (file.exists("D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/tm/sea_ice_data.csv")) {
  sample_output <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/tm/sea_ice_data.csv")
  print(head(sample_output))
}
