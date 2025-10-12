# Load necessary libraries
library(terra)
library(sf)
library(dplyr)
library(stringr)

# Function to compute mean duration, standard deviation, mean persistence, and persistence standard deviation of sea ice concentration above a threshold
compute_duration_persistence_stats <- function(buffer_path, nsidc, threshold = .15, winter_months = c(6, 7, 8, 9)) {
  combined_buffers_sf <- st_read(buffer_path)
  dissolved_buffer <- st_union(combined_buffers_sf)
  
  # Extract dates before masking
  dates <- time(nsidc)
  
  buffer_mask <- mask(nsidc, vect(dissolved_buffer))
  
  all_years <- unique(year(dates))
  
  duration_persistence_stats <- data.frame(
    year = integer(),
    month = character(),
    mean_duration = numeric(),
    sd_duration = numeric(),
    mean_persistence = numeric(),
    sd_persistence = numeric()
  )
  
  calculate_metrics <- function(data, threshold) {
    durations <- app(data, function(x) {
      rle_result <- rle(x > threshold)
      durations <- rle_result$lengths[rle_result$values]
      return(mean(durations, na.rm = TRUE))
    })
    mean_duration <- mean(values(durations), na.rm = TRUE)
    sd_duration <- sd(values(durations), na.rm = TRUE)
    
    open_water_prop <- app(data, function(x) mean(x < threshold, na.rm = TRUE))
    mean_persistence <- mean(values(open_water_prop), na.rm = TRUE)
    sd_persistence <- sd(values(open_water_prop), na.rm = TRUE)
    
    return(list(mean_duration = mean_duration, sd_duration = sd_duration, mean_persistence = mean_persistence, sd_persistence = sd_persistence))
  }
  
  for (year in all_years) {
    for (month in winter_months) {
      monthly_indices <- which(year(dates) == year & month(dates) == month)
      if (length(monthly_indices) > 0) {
        monthly_data <- buffer_mask[[monthly_indices]]
        metrics <- calculate_metrics(monthly_data, threshold)
        
        duration_persistence_stats <- rbind(duration_persistence_stats, data.frame(
          year = year,
          month = month,
          mean_duration = metrics$mean_duration,
          sd_duration = metrics$sd_duration,
          mean_persistence = metrics$mean_persistence,
          sd_persistence = metrics$sd_persistence
        ))
      }
    }
    
    season_indices <- which(year(dates) == year & month(dates) %in% winter_months)
    if (length(season_indices) > 0) {
      season_data <- buffer_mask[[season_indices]]
      metrics <- calculate_metrics(season_data, threshold)
      
      duration_persistence_stats <- rbind(duration_persistence_stats, data.frame(
        year = year,
        month = "Season-wide",
        mean_duration = metrics$mean_duration,
        sd_duration = metrics$sd_duration,
        mean_persistence = metrics$mean_persistence,
        sd_persistence = metrics$sd_persistence
      ))
    }
  }
  
  return(duration_persistence_stats)
}

# Main analysis function
analyze_sea_ice_effect <- function(thresholds, nsidc_subset, results_dir) {
  # Initialize dataframe to store all results
  all_duration_persistence_stats <- data.frame()
  
  # Define home range directory
  home_range_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/home-ranges"
  home_range_files <- list.files(home_range_dir, pattern = "\\.shp$", full.names = TRUE)
  
  # Loop through each home range shapefile and compute results
  for (threshold in thresholds) {
    for (buffer_path in home_range_files) {
      home_range_size <- str_extract(buffer_path, "(\\d+)km")
      duration_persistence_stats <- compute_duration_persistence_stats(buffer_path, nsidc_subset, threshold = threshold)
      duration_persistence_stats$Threshold <- threshold
      duration_persistence_stats$HomeRangeSize <- home_range_size
      all_duration_persistence_stats <- bind_rows(all_duration_persistence_stats, duration_persistence_stats)
    }
  }
  
  # Export the compiled results to CSV
  write.csv(all_duration_persistence_stats, file.path(results_dir, "sea_ice_duration_persistence_stats_subset.csv"), row.names = FALSE)
  
  # Export the compiled results to RDS
  saveRDS(all_duration_persistence_stats, file.path(results_dir, "sea_ice_duration_persistence_stats.rds"))
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
analyze_sea_ice_effect(c(.15, .30, .50), nsidc, results_dir)

# Display the results
all_duration_persistence_stats <- read.csv(file.path(results_dir, "sea_ice_duration_persistence_stats.csv"))
head(all_duration_persistence_stats)
