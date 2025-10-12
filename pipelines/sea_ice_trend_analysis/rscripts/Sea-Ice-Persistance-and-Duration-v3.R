# Load necessary libraries
library(terra)
library(dplyr)
library(lubridate)
library(nlme)
library(ggplot2)
library(grDevices)  # For colorRampPalette

# Define file paths
chunk_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/25km_Sea-Ice-Index/stack/substack"
chunk_files <- list.files(chunk_dir, pattern = "\\.nc$", full.names = TRUE)
study_area_shapefile <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/shp/subregions/Frostbound_AQ_Subregions_EPSG_3976.shp"

# Load the study area shapefile
study_area <- vect(study_area_shapefile)

# Define the threshold for sea ice concentration
ice_threshold <- 15

# Define the date range for analysis
start_year <- 1981
end_year <- 2023

# Function to calculate consecutive days above threshold
calculate_duration <- function(x, threshold) {
  if (all(is.na(x))) return(NA)
  rle_result <- rle(x > threshold)
  max_duration <- ifelse(any(rle_result$values), max(rle_result$lengths[rle_result$values]), 0)
  return(max_duration)
}

# Simplified function to process each chunk file using app
process_chunk <- function(file_path) {
  region_name <- gsub(".*/|\\.nc$", "", file_path)
  message(paste("Processing file:", file_path))
  
  # Load the NSIDC sea ice concentration data for the subregion
  nsidc <- rast(file_path)
  
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

# Mask the rasters by the study area
persistence_stack_masked <- mask(persistence_stack, study_area)
duration_stack_masked <- mask(duration_stack, study_area)

# Function to apply linear model to each pixel and return slope and p-value
apply_lm <- function(y) {
  if (all(is.na(y))) {
    return(c(NA, NA))
  } else {
    years <- start_year:end_year
    fit <- lm(y ~ years)
    slope <- coef(fit)[2]  # Slope of the linear model
    p_value <- summary(fit)$coefficients[2, 4]  # p-value for the slope
    return(c(slope, p_value))
  }
}

# Apply the linear model to each pixel in the masked persistence stack
persistence_trend <- app(persistence_stack_masked, apply_lm)

# Apply the linear model to each pixel in the masked duration stack
duration_trend <- app(duration_stack_masked, apply_lm)

# Extract slopes and p-values from the results
persistence_slope <- persistence_trend[[1]]
persistence_pvalue <- persistence_trend[[2]]
duration_slope <- duration_trend[[1]]
duration_pvalue <- duration_trend[[2]]

# Verify the integrity of the data
print("Summary of persistence slope layer:")
print(summary(persistence_slope))

print("Summary of persistence p-value layer:")
print(summary(persistence_pvalue))

print("Summary of duration slope layer:")
print(summary(duration_slope))

print("Summary of duration p-value layer:")
print(summary(duration_pvalue))

# Define a significance level
alpha <- 0.05

# Calculate statistics for persistence
persistence_values <- values(persistence_slope)
persistence_pvalues <- values(persistence_pvalue)

# Remove NA values
persistence_data <- cbind(persistence_values, persistence_pvalues)
persistence_data <- persistence_data[complete.cases(persistence_data), ]
persistence_values <- persistence_data[, 1]
persistence_pvalues <- persistence_data[, 2]

# Calculate mean, median, min, max slopes for persistence
persistence_mean <- mean(persistence_values)
persistence_median <- median(persistence_values)
persistence_min <- min(persistence_values)
persistence_max <- max(persistence_values)

# Calculate percentage of significant positive, significant negative, and non-significant trends
persistence_significant <- sum(persistence_pvalues <= alpha) / length(persistence_pvalues) * 100
persistence_positive <- sum(persistence_pvalues <= alpha & persistence_values > 0) / length(persistence_pvalues) * 100
persistence_negative <- sum(persistence_pvalues <= alpha & persistence_values < 0) / length(persistence_pvalues) * 100
persistence_nonsignificant <- 100 - persistence_significant

# Calculate statistics for duration
duration_values <- values(duration_slope)
duration_pvalues <- values(duration_pvalue)

# Remove NA values
duration_data <- cbind(duration_values, duration_pvalues)
duration_data <- duration_data[complete.cases(duration_data), ]
duration_values <- duration_data[, 1]
duration_pvalues <- duration_data[, 2]

# Calculate mean, median, min, max slopes for duration
duration_mean <- mean(duration_values)
duration_median <- median(duration_values)
duration_min <- min(duration_values)
duration_max <- max(duration_values)

# Calculate percentage of significant positive, significant negative, and non-significant trends
duration_significant <- sum(duration_pvalues <= alpha) / length(duration_pvalues) * 100
duration_positive <- sum(duration_pvalues <= alpha & duration_values > 0) / length(duration_pvalues) * 100
duration_negative <- sum(duration_pvalues <= alpha & duration_values < 0) / length(duration_pvalues) * 100
duration_nonsignificant <- 100 - duration_significant

# Print summary statistics for reporting
persistence_summary <- data.frame(
  Metric = c("Mean slope", "Median slope", "Minimum slope", "Maximum slope", "% Significant", "% Positive", "% Negative", "% Non-significant"),
  Persistence = c(persistence_mean, persistence_median, persistence_min, persistence_max, persistence_significant, persistence_positive, persistence_negative, persistence_nonsignificant)
)

duration_summary <- data.frame(
  Metric = c("Mean slope", "Median slope", "Minimum slope", "Maximum slope", "% Significant", "% Positive", "% Negative", "% Non-significant"),
  Duration = c(duration_mean, duration_median, duration_min, duration_max, duration_significant, duration_positive, duration_negative, duration_nonsignificant)
)

print(persistence_summary)
print(duration_summary)

# Define a color palette for trends
trend_colors <- colorRampPalette(c("red", "white", "blue"))

# Classify the trends based on significance and directionality using ifel function
persistence_significance_direction <- ifel(persistence_pvalue <= alpha & persistence_slope > 0, 1,
                                           ifel(persistence_pvalue <= alpha & persistence_slope < 0, -1, 0))

duration_significance_direction <- ifel(duration_pvalue <= alpha & duration_slope > 0, 1,
                                        ifel(duration_pvalue <= alpha & duration_slope < 0, -1, 0))

# Define a color palette for significance and directionality
significance_colors <- c("red", "white", "blue")  # Red for negative significant, white for non-significant, blue for positive significant

# Plotting the trends with custom colors
plot(persistence_slope, main="Trend of Persistence Over Time", col=trend_colors(100))
plot(duration_slope, main="Trend of Duration Over Time", col=trend_colors(100))

# Plotting significance and directionality
plot(persistence_significance_direction, main="Significance and Directionality of Persistence Trends", col=significance_colors, legend=FALSE)
plot(duration_significance_direction, main="Significance and Directionality of Duration Trends", col=significance_colors, legend=FALSE)

# Additional Verification of Summary Metrics
print("Verification of negative trends in persistence and duration:")
print(sum(persistence_values < 0 & persistence_pvalues <= alpha))
print(sum(duration_values < 0 & duration_pvalues <= alpha))
