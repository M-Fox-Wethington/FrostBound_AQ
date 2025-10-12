# Load necessary libraries
library(terra)
library(dplyr)
library(lubridate)
library(nlme)

# Define the directory containing the subregion files
chunk_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/25km_Sea-Ice-Index/stack/substack"
chunk_files <- list.files(chunk_dir, pattern = "\\.nc$", full.names = TRUE)

# Define the threshold for sea ice concentration
ice_threshold <- 15

# Define the date range for analysis
start_year <- 1981
end_date <- as.Date("2023-09-30")

# Function to calculate consecutive days above threshold
calculate_duration <- function(x, threshold) {
  if (all(is.na(x))) return(NA)  # Return NA if all values are NA
  rle_result <- rle(x > threshold)
  max_duration <- ifelse(any(rle_result$values), max(rle_result$lengths[rle_result$values]), 0)
  return(max_duration)
}

# Simplified function to process each chunk file using app
process_chunk <- function(file_path) {
  region_name <- gsub(".*/|\\.nc$", "", file_path)  # Extract the region name from the file path
  message(paste("Processing file:", file_path))
  
  # Load the NSIDC sea ice concentration data for the subregion
  nsidc <- rast(file_path)
  print("Loaded NSIDC data:")
  print(nsidc)
  
  # Filter the sea ice data from 1981 to 2023-09-30
  nsidc <- subset(nsidc, which(time(nsidc) >= as.Date(paste0(start_year, "-01-01")) & time(nsidc) <= end_date))
  print("Filtered NSIDC data:")
  print(nsidc)
  
  # Calculate persistence: number of days with ice concentration above the threshold for each pixel
  persistence <- app(nsidc, function(x) sum(x > ice_threshold, na.rm = TRUE))
  names(persistence) <- "Persistence"
  
  # Calculate duration: maximum consecutive days with ice concentration above the threshold for each pixel
  duration <- app(nsidc, function(x) calculate_duration(x, ice_threshold))
  names(duration) <- "Duration"
  
  # Calculate the statistics for each annual raster layer
  persistence_stats <- data.frame(
    Region = region_name,
    Min_Persistence = min(persistence[], na.rm = TRUE),
    Max_Persistence = max(persistence[], na.rm = TRUE),
    Mean_Persistence = mean(persistence[], na.rm = TRUE),
    Median_Persistence = median(persistence[], na.rm = TRUE),
    Start_Year = start_year,
    End_Date = end_date
  )
  
  duration_stats <- data.frame(
    Region = region_name,
    Min_Duration = min(duration[], na.rm = TRUE),
    Max_Duration = max(duration[], na.rm = TRUE),
    Mean_Duration = mean(duration[], na.rm = TRUE),
    Median_Duration = median(duration[], na.rm = TRUE),
    Start_Year = start_year,
    End_Date = end_date
  )
  
  # Plot the rasters
  plot(persistence, main = paste("Sea Ice Persistence Above Threshold for", region_name), col = hcl.colors(100, "Blues", rev = TRUE))
  plot(duration, main = paste("Sea Ice Duration Above Threshold for", region_name), col = hcl.colors(100, "Reds", rev = TRUE))
  
  return(list(persistence_stats, duration_stats))
}

# Process each chunk file sequentially
results <- lapply(chunk_files, process_chunk)

# Combine all results into single dataframes
persistence_df <- do.call(rbind, lapply(results, `[[`, 1))
duration_df <- do.call(rbind, lapply(results, `[[`, 2))

# View the results
print(persistence_df)
print(duration_df)

# Fit GLS models to the calculated statistics
perform_gls <- function(stats_df, metric) {
  gls_model <- gls(as.formula(paste(metric, "~ Start_Year")), data = stats_df, correlation = corAR1())
  summary(gls_model)
}

# Fit GLS models for each region and metric
gls_results_persistence <- lapply(unique(persistence_df$Region), function(region) {
  region_stats <- subset(persistence_df, Region == region)
  gls_result <- perform_gls(region_stats, "Mean_Persistence")
  return(gls_result)
})

gls_results_duration <- lapply(unique(duration_df$Region), function(region) {
  region_stats <- subset(duration_df, Region == region)
  gls_result <- perform_gls(region_stats, "Mean_Duration")
  return(gls_result)
})

# Print GLS summaries
print(gls_results_persistence)
print(gls_results_duration)

# Save the results to files
write.csv(persistence_df, "persistence_df.csv")
write.csv(duration_df, "duration_df.csv")
