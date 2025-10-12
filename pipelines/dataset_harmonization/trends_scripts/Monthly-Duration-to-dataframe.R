# Load required libraries
library(terra)
library(dplyr)
library(lubridate)

# Define the function to calculate duration for a given raster stack
calculate_duration_stats <- function(file_path, ice_threshold = 0.15, winter_months = c(6, 7, 8, 9)) {
  # Load the raster stack from file
  raster_stack <- rast(file_path)
  
  # Extract dates from the raster stack
  layer_dates <- as.Date(as.numeric(time(raster_stack)), origin = "1970-01-01")
  
  # Filter for winter months (June - September)
  winter_indices <- which(month(layer_dates) %in% winter_months)
  winter_raster_stack <- subset(raster_stack, winter_indices)
  winter_dates <- layer_dates[winter_indices]
  
  # Calculate duration: maximum length of consecutive days with SIC >= threshold for each month
  calculate_duration <- function(x, threshold) {
    rle_result <- rle(as.vector(x) >= threshold)
    max_duration <- ifelse(any(rle_result$values), max(rle_result$lengths[rle_result$values]), 0)
    return(max_duration)
  }
  
  months <- unique(floor_date(winter_dates, "month"))
  duration_list <- vector("list", length(months))
  names(duration_list) <- months
  
  for (month in months) {
    month_indices <- which(floor_date(winter_dates, "month") == month)
    month_raster_stack <- subset(winter_raster_stack, month_indices)
    duration <- app(month_raster_stack, function(x) calculate_duration(x, ice_threshold))
    duration_list[[as.character(month)]] <- global(duration, fun = mean, na.rm = TRUE)[, 1]
  }
  
  # Extract region name from the file path
  region_name <- tools::file_path_sans_ext(basename(file_path))
  region_name <- sub("NSIDC_25km_Harmonized_", "", region_name)
  
  # Combine results into a dataframe
  duration_df <- data.frame(
    Month = as.Date(names(duration_list)),
    Mean_Monthly_Duration = unlist(duration_list),
    Region = region_name
  )
  
  return(duration_df)
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
  
  # Run the calculate_duration_stats function on the current file
  result <- calculate_duration_stats(file)
  
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
# write.csv(master_df, file = "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/trend-analysis/Duration_Metrics.csv", row.names = FALSE)

# Save the master dataframe as RDS (optional)
# saveRDS(master_df, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/trend-analysis/Duration_Metrics.rds")
