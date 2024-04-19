library(terra)

# Define the file path for the raster stack
r_fp <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/AMSR-Unified_SIC_Full-Catalog_Time_Series.nc"

# Load the raster stack outside of the loop
raster_stack <- rast(r_fp)

# Function to calculate mean and check for sea ice extent
calculate_layer_stats <- function(layer_data) {
  # Cast values of 0 to NA
  values(layer_data)[values(layer_data) == 0] <- NA
  
  # Calculate mean and sea ice extent only if not all values are NA
  mean_value <- NA
  ice_extent <- NA
  if (any(!is.na(values(layer_data)))) {
    mean_value <- mean(values(layer_data), na.rm = TRUE)
    # Convert sea ice concentration to binary and calculate extent
    binary_values <- ifelse(values(layer_data) >= 0.15, 1, 0)
    ice_extent <- sum(binary_values * res(layer_data)[1] * res(layer_data)[2], na.rm = TRUE)
  }
  return(list(mean = mean_value, ice_extent = ice_extent))
}

# Extract dates for each layer
dates <- as.Date(time(raster_stack))

# Initialize result list
results <- vector("list", length = nlyr(raster_stack))
names(results) <- 1:nlyr(raster_stack)

# Perform calculations sequentially
for (i in 1:nlyr(raster_stack)) {
  print(paste("Processing Layer:", i))
  layer_data <- raster_stack[[i]]  # Extract layer from stack
  results[[i]] <- calculate_layer_stats(layer_data)
  rm(layer_data)  # Remove the temporary layer data
  gc()  # Explicitly call garbage collection
}

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
