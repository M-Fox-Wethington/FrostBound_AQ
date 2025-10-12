
# This script analyzes daily sea ice concentration (SIC) and extent within defined home ranges by calculating the daily mean SIC, standard deviation of SIC, and total ice-covered area (sea ice extent) above a threshold for each day. It loads NSIDC sea ice data, applies a threshold to set low values to zero, and computes the specified metrics for each home range. The sea ice extent is calculated as the total area of cells with SIC above the threshold, effectively representing the "total ice-covered area" within the home range. The results are compiled and exported to a CSV file for further analysis.

# Load necessary libraries
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(lubridate)

# Function to compute daily mean SIC, SD SIC, and sea ice extent above a threshold
compute_daily_sic_statistics <- function(buffer_path, nsidc, threshold = 0.15, cell_area_sq_meters) {
  combined_buffers_sf <- st_read(buffer_path)
  dissolved_buffer <- st_union(combined_buffers_sf)
  
  # Extract dates before masking
  dates <- time(nsidc)
  
  buffer_mask <- mask(nsidc, vect(dissolved_buffer))
  
  # Set all values < threshold to 0 for the entire raster stack
  buffer_mask <- app(buffer_mask, fun = function(x) { x[x < threshold] <- 0; return(x) })
  
  # Calculate mean and SD SIC excluding NAs for each layer
  mean_sic <- global(buffer_mask, fun = 'mean', na.rm = TRUE)[, 1]
  sd_sic <- global(buffer_mask, fun = 'sd', na.rm = TRUE)[, 1]
  
  # Calculate sea ice extent (total cells above threshold in square kilometers)
  valid_ice_cells <- global(buffer_mask >= threshold, fun = 'sum', na.rm = TRUE)[, 1]
  total_ice_area_sq_km <- (valid_ice_cells * cell_area_sq_meters) / 1e6  # Convert total ice area to square kilometers
  
  results <- data.frame(
    date = as.Date(dates, origin = "1970-01-01"),
    mean_sic = mean_sic,
    sd_sic = sd_sic,
    ice_extent_km2 = total_ice_area_sq_km
  )
  
  return(results)
}

# Main analysis function
analyze_sea_ice_effect <- function(thresholds, nsidc_subset, results_dir) {
  # Initialize dataframe to store all results
  all_results_df <- data.frame()
  
  # Calculate the cell area in square meters (only once)
  cell_area_sq_meters <- prod(res(nsidc_subset))
  
  # Loop through each home range shapefile and compute results
  home_range_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/home-ranges"
  home_range_files <- list.files(home_range_dir, pattern = "\\.shp$", full.names = TRUE)
  
  for (threshold in thresholds) {
    for (buffer_path in home_range_files) {
      home_range_size <- str_extract(buffer_path, "(\\d+)km")
      metrics <- compute_daily_sic_statistics(buffer_path, nsidc_subset, threshold = threshold, cell_area_sq_meters = cell_area_sq_meters)
      
      metrics <- metrics %>%
        mutate(Threshold = threshold, HomeRangeSize = home_range_size)
      
      all_results_df <- bind_rows(all_results_df, metrics)
    }
  }
  
  # Export the compiled results to CSV
  write.csv(all_results_df, file.path(results_dir, "daily_sic_statistics.csv"), row.names = FALSE)
  
  return(all_results_df)
}

# Load the NSIDC sea ice concentration data
nsidc <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/complete-harmonized-dataset/tif/nsidc_12_5km_harmonized_1979-2023.tif")

# Filter the sea ice data from 1981 to 2023
start_date <- as.Date("1981-01-01")
end_date <- as.Date("2023-09-30")
nsidc <- subset(nsidc, which(time(nsidc) >= start_date & time(nsidc) <= end_date))

# Subset the NSIDC data for testing (e.g., first 100 layers)
# nsidc_subset <- subset(nsidc, 1:100)

# Create results directory
results_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/metric-calculation-csv"
dir.create(results_dir, showWarnings = FALSE)

# Run the analysis on the subset
all_results_df <- analyze_sea_ice_effect(c(0.15), nsidc, results_dir)

# Display the results
head(all_results_df)


