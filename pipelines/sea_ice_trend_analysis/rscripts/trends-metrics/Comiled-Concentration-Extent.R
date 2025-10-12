# Load necessary libraries
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(lubridate)

# Function to compute extent statistics of sea ice concentration above a threshold
compute_extent_statistics <- function(buffer_path, nsidc, threshold = 0.15, winter_months = c(6, 7, 8, 9)) {
  combined_buffers_sf <- st_read(buffer_path)
  dissolved_buffer <- st_union(combined_buffers_sf)
  buffer_mask <- mask(nsidc, vect(dissolved_buffer))
  all_years <- unique(year(time(buffer_mask)))
  
  results <- data.frame(
    year = rep(all_years, each = 6),
    period = rep(c("annual", "winter", "June", "July", "August", "September"), length(all_years)),
    min_extent = NA,
    max_extent = NA,
    mean_extent = NA
  )
  
  cell_area_sq_meters <- prod(res(buffer_mask[[1]]))
  
  for (year in all_years) {
    annual_data <- subset(buffer_mask, which(year(time(buffer_mask)) == year))
    winter_data <- subset(buffer_mask, which(year(time(buffer_mask)) == year & month(time(buffer_mask)) %in% winter_months))
    monthly_data <- lapply(winter_months, function(m) subset(buffer_mask, which(year(time(buffer_mask)) == year & month(time(buffer_mask)) == m)))
    
    periods_data <- list(annual = annual_data, winter = winter_data, June = monthly_data[[1]], July = monthly_data[[2]], August = monthly_data[[3]], September = monthly_data[[4]])
    
    for (period in names(periods_data)) {
      data <- periods_data[[period]]
      if (nlyr(data) > 0) {
        # Set all values < threshold to 0
        data <- app(data, fun = function(x) { x[x < threshold] <- 0; return(x) })
        
        # Calculate extent excluding NAs for each layer
        mean_extent <- global(data, fun = 'mean', na.rm = TRUE)[, 1]
        min_extent <- global(data, fun = 'min', na.rm = TRUE)[, 1]
        max_extent <- global(data, fun = 'max', na.rm = TRUE)[, 1]
        
        row_index <- which(results$year == year & results$period == period)
        results$mean_extent[row_index] <- mean(mean_extent, na.rm = TRUE)
        results$min_extent[row_index] <- min(min_extent, na.rm = TRUE)
        results$max_extent[row_index] <- max(max_extent, na.rm = TRUE)
      }
    }
  }
  
  return(results)
}

# Main analysis function
analyze_sea_ice_effect <- function(thresholds) {
  # Load Gentoo penguin data
  penguin_data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/inputs/gentoo_presence_absence_assumptions.csv")
  penguin_abundance_data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/inputs/modeled_gentoo_parameters.csv")
  
  # Load study area shapefile
  study_area_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/shp/Frostbound_Study_Areas_EPSG_3976.shp"
  study_area <- st_read(study_area_path)
  
  # Load the NSIDC sea ice concentration data
  nsidc <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/complete-harmonized-dataset/tif/nsidc_12_5km_harmonized_1979-2023.tif")
  
  # Filter the sea ice data from 1981 to 2023
  start_date <- as.Date("1981-01-01")
  end_date <- as.Date("2023-09-30")
  nsidc <- subset(nsidc, which(time(nsidc) >= start_date & time(nsidc) <= end_date))
  
  # Create results directory
  results_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/calculated-metrics"
  dir.create(results_dir, showWarnings = FALSE)
  
  # Initialize dataframe to store all results
  all_results_df <- data.frame()
  
  # Loop through each home range shapefile and compute results
  home_range_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/home-ranges"
  home_range_files <- list.files(home_range_dir, pattern = "\\.shp$", full.names = TRUE)
  
  for (threshold in thresholds) {
    for (buffer_path in home_range_files) {
      home_range_size <- str_extract(buffer_path, "(\\d+)km")
      metrics <- compute_extent_statistics(buffer_path, nsidc, threshold = threshold)
      
      metrics <- metrics %>%
        mutate(Threshold_Extent = threshold, HomeRangeSize_Extent = home_range_size)
      
      all_results_df <- rbind(all_results_df, metrics)
    }
  }
  
  # Export the compiled results to CSV
  write.csv(all_results_df, file.path(results_dir, "sea_ice_extent_statistics.csv"), row.names = FALSE)
  
  # Optionally, save to RDS
  saveRDS(all_results_df, file.path(results_dir, "sea_ice_extent_statistics.rds"))
}

# Run the analysis
analyze_sea_ice_effect(c(0.15))
