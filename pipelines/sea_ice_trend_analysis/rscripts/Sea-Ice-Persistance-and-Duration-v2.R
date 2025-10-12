# Load necessary libraries
library(terra)
library(dplyr)
library(lubridate)
library(nlme)
library(ggplot2)

# Define the directory containing the subregion files
chunk_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/25km_Sea-Ice-Index/stack/substack"
chunk_files <- list.files(chunk_dir, pattern = "\\.nc$", full.names = TRUE)

# Define the threshold for sea ice concentration
ice_threshold <- 15

# Define the date range for analysis
start_year <- 1981
end_year <- 2023

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
  
  annual_persistence <- list()
  annual_duration <- list()
  
  for (year in start_year:end_year) {
    message(paste("Processing year:", year))
    start_date <- as.Date(paste0(year, "-01-01"))
    end_date <- as.Date(paste0(year, "-12-31"))
    
    # Filter the sea ice data for the current year
    nsidc_year <- subset(nsidc, which(time(nsidc) >= start_date & time(nsidc) <= end_date))
    
    if (nlyr(nsidc_year) == 0) {
      message(paste("No data for year:", year))
      next
    }
    
    # Calculate persistence: number of days with ice concentration above the threshold for each pixel
    persistence <- app(nsidc_year, function(x) sum(x > ice_threshold, na.rm = TRUE))
    names(persistence) <- paste0("Persistence_", year)
    message("Calculated persistence")
    
    # Calculate duration: maximum consecutive days with ice concentration above the threshold for each pixel
    duration <- app(nsidc_year, function(x) calculate_duration(x, ice_threshold))
    names(duration) <- paste0("Duration_", year)
    message("Calculated duration")
    
    annual_persistence[[as.character(year)]] <- persistence
    annual_duration[[as.character(year)]] <- duration
  }
  
  message(paste("Finished processing file:", file_path))
  return(list(persistence = annual_persistence, duration = annual_duration))
}

# Process each chunk file sequentially
results <- lapply(chunk_files, process_chunk)

# Combine all annual rasters into a single stack for persistence and duration
combine_rasters <- function(raster_list) {
  raster_stack <- rast()
  for (year in names(raster_list)) {
    if (!is.null(raster_list[[year]])) {
      raster_stack <- c(raster_stack, raster_list[[year]])
    }
  }
  return(raster_stack)
}

# Combine results for all files
persistence_stack <- combine_rasters(do.call(c, lapply(results, `[[`, "persistence")))
duration_stack <- combine_rasters(do.call(c, lapply(results, `[[`, "duration")))

# Check if stacks are created correctly
print(persistence_stack)
print(duration_stack)

# Function to fit linear model to each pixel and extract the slope
calculate_trend <- function(raster_stack) {
  years <- as.numeric(sub(".*_(\\d{4})$", "\\1", names(raster_stack)))
  trend_raster <- rast(raster_stack[[1]])
  values(trend_raster) <- NA
  
  for (row in 1:nrow(raster_stack)) {
    for (col in 1:ncol(raster_stack)) {
      pixel_values <- as.numeric(values(raster_stack, row = row, col = col))
      if (all(is.na(pixel_values))) next
      
      lm_result <- try(lm(pixel_values ~ years), silent = TRUE)
      if (inherits(lm_result, "try-error")) next
      
      slope <- coef(lm_result)[2]
      trend_raster[row, col] <- slope
    }
  }
  return(trend_raster)
}

# Calculate trends for persistence and duration
persistence_trend <- calculate_trend(persistence_stack)
duration_trend <- calculate_trend(duration_stack)

# Plot the trends
plot(persistence_trend, main = "Trend in Sea Ice Persistence", col = hcl.colors(100, "Blues", rev = TRUE))
plot(duration_trend, main = "Trend in Sea Ice Duration", col = hcl.colors(100, "Reds", rev = TRUE))

# Save the results to files
writeRaster(persistence_trend, "persistence_trend.tif", overwrite = TRUE)
writeRaster(duration_trend, "duration_trend.tif", overwrite = TRUE)
