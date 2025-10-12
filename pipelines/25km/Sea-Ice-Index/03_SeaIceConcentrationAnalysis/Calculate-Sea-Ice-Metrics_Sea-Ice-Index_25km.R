# This script processes chunks of satellite-derived sea ice concentration data to calculate 
# mean Sea Ice Concentration (SIC) and total sea ice extent for each raster layer within the chunks.
# Each chunk file represents a portion of a comprehensive sea ice index catalog, starting from 1978 to the present.
# The script utilizes parallel computing to efficiently handle computations for each layer within the chunk files.
# Results are compiled into a CSV and RDS file for further analysis or reporting.

library(terra)        # For raster data manipulation
library(doParallel)   # For parallel processing
library(foreach)      # For looping with parallel support

# Function to calculate mean SIC and sea ice extent for a given raster layer within a chunk file
calculate_layer_stats <- function(layer_index, file_path) {
  raster_stack <- rast(file_path)                  # Load the raster stack from file
  layer_data <- raster_stack[[layer_index]]        # Extract the specific layer
  
  layer_date <- time(layer_data)                   # Extract date from layer metadata
  
  mean_sic <- as.numeric(terra::global(layer_data, fun = 'mean', na.rm = TRUE))  # Calculate mean SIC excluding NAs
  
  cell_area_sq_meters <- prod(res(layer_data))    # Calculate the area of one cell in square meters
  
  valid_ice_cells <- sum(values(layer_data) >= 15.0, na.rm = TRUE)  # Count cells with valid ice concentration
  
  total_ice_area_sq_km <- (valid_ice_cells * cell_area_sq_meters) / 1e6  # Convert total ice area to square kilometers
  
  rm(raster_stack, layer_data)   # Clear memory
  gc()                           # Run garbage collection
  
  return(list(mean_sic = mean_sic, ice_extent_km = total_ice_area_sq_km, date = layer_date))
}

chunk_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/raster-stack"  # Directory with chunk files
chunk_files <- list.files(chunk_dir, pattern = "\\.nc$", full.names = TRUE)  # List all .nc files

all_results <- list()  # Initialize list to store results from all chunks

for (r_fp in chunk_files) {
  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("calculate_layer_stats", "r_fp"))
  
  results <- foreach(i = 1:nlyr(rast(r_fp)), .packages = 'terra', .errorhandling = "pass") %dopar% {
    tryCatch({
      calculate_layer_stats(i, r_fp)
    }, error = function(e) {
      return(list(error = TRUE, message = e$message))
    })
  }
  
  stopCluster(cl)
  
  all_results <- c(all_results, results)  # Append chunk results to overall results
}

# Convert list of all results into a dataframe, handling possible errors
final_result_df <- do.call(rbind, lapply(all_results, function(x) {
  if (!is.null(x$error)) {
    return(data.frame(Date = NA, MeanSIC = NA, IceExtent_km = NA, error = x$message))
  } else {
    return(data.frame(Date = x$date, MeanSIC = x$mean_sic, IceExtent_km = x$ice_extent_km, error = NA))
  }
}))

# Save combined results to CSV and RDS for further use
write.csv(final_result_df, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/tm/sea_ice_final_data.csv", row.names = FALSE)
saveRDS(final_result_df, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/tm/sea_ice_final_data.rds")
