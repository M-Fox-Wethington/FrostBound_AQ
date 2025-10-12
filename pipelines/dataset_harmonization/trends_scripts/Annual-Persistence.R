# Load required libraries
library(terra)
library(dplyr)
library(lubridate)

# Define the function to calculate persistence for a given raster stack
calculate_persistence_stats <- function(file_path, ice_threshold = 0.15, winter_months = c(6, 7, 8, 9)) {
  # Load the raster stack from file
  raster_stack <- rast(file_path)
  
  # Extract dates from the raster stack
  layer_dates <- as.Date(as.numeric(time(raster_stack)), origin = "1970-01-01")
  
  # Filter for winter months (June - September)
  winter_indices <- which(month(layer_dates) %in% winter_months)
  winter_raster_stack <- subset(raster_stack, winter_indices)
  winter_dates <- layer_dates[winter_indices]
  
  # Calculate persistence: sum of days with SIC >= threshold for each year
  calculate_persistence <- function(x, threshold) {
    rle_result <- rle(as.vector(x) >= threshold)
    max_duration <- ifelse(any(rle_result$values), max(rle_result$lengths[rle_result$values]), 0)
    return(max_duration)
  }
  
  years <- unique(year(winter_dates))
  persistence_list <- vector("list", length(years))
  names(persistence_list) <- years
  
  for (year in years) {
    year_indices <- which(year(winter_dates) == year)
    year_raster_stack <- subset(winter_raster_stack, year_indices)
    persistence <- app(year_raster_stack, function(x) calculate_persistence(x, ice_threshold))
    persistence_list[[as.character(year)]] <- global(persistence, fun = mean, na.rm = TRUE)[, 1]
  }
  
  # Extract region name from the file path
  region_name <- tools::file_path_sans_ext(basename(file_path))
  region_name <- sub("NSIDC_25km_Harmonized_", "", region_name)
  
  # Combine results into a dataframe
  persistence_df <- data.frame(
    Year = as.integer(names(persistence_list)),
    Persistence = unlist(persistence_list),
    Region = region_name
  )
  
  return(persistence_df)
}

# Directory containing the .tif files
chunk_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/complete-harmonized-dataset/tif"

# List all .tif files for each region
chunk_files <- list.files(chunk_dir, pattern = "\\.tif$", full.names = TRUE)

# Initialize an empty dataframe to store the combined results
master_df <- data.frame()

# Iterate through each .tif file and combine the results into master_df
for (file in chunk_files) {
  cat("Processing file:", file, "\n")
  
  # Run the calculate_persistence_stats function on the current file
  result <- calculate_persistence_stats(file)
  
  # Combine the result into the master dataframe
  master_df <- bind_rows(master_df, result)
}

# Rename the region "nsidc_12_5km_harmonized_1979-2023" to "Complete Study Area"
master_df <- master_df %>%
  mutate(Region = ifelse(Region == "nsidc_12_5km_harmonized_1979-2023", "Complete Study Area", Region))

# Print the structure of the master dataframe
str(master_df)

# Print the first few rows of the master dataframe
head(master_df)

# Save the master dataframe to a CSV file (optional)
# write.csv(master_df, file = "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/trend-analysis/Persistence_Metrics.csv", row.names = FALSE)

# Save the master dataframe as RDS (optional)
# saveRDS(master_df, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/trend-analysis/Persistence_Metrics.rds")
