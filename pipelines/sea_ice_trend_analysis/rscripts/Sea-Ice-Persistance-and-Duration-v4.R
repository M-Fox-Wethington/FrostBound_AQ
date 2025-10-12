# Load necessary libraries
library(terra)
library(dplyr)
library(lubridate)
library(nlme)
library(ggplot2)
library(grDevices)  # For colorRampPalette
library(gridExtra)  # For arranging plots side by side

# Define the function
analyze_sea_ice <- function(start_year, end_year, use_mask = TRUE, report_metrics = c("concentration", "duration", "persistence")) {
  
  # Define file paths
  chunk_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/25km_Sea-Ice-Index/stack/substack"
  chunk_files <- list.files(chunk_dir, pattern = "\\.nc$", full.names = TRUE)
  study_area_shapefile <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/shp/subregions/Frostbound_AQ_Subregions_EPSG_3976.shp"
  
  # Load the study area shapefile
  if (use_mask) {
    study_area <- vect(study_area_shapefile)
  }
  
  # Define the threshold for sea ice concentration
  ice_threshold <- 15
  
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
    annual_concentration <- list()
    
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
      
      # Calculate metrics based on options
      if ("persistence" %in% report_metrics) {
        persistence <- app(nsidc_year, function(x) sum(x > ice_threshold, na.rm = TRUE))
        names(persistence) <- paste0("Persistence_", year)
        annual_persistence[[as.character(year)]] <- persistence
      }
      
      if ("duration" %in% report_metrics) {
        duration <- app(nsidc_year, function(x) calculate_duration(x, ice_threshold))
        names(duration) <- paste0("Duration_", year)
        annual_duration[[as.character(year)]] <- duration
      }
      
      if ("concentration" %in% report_metrics) {
        concentration <- app(nsidc_year, function(x) mean(x, na.rm = TRUE))
        names(concentration) <- paste0("Concentration_", year)
        annual_concentration[[as.character(year)]] <- concentration
      }
    }
    
    message(paste("Finished processing file:", file_path))
    return(list(persistence = annual_persistence, duration = annual_duration, concentration = annual_concentration))
  }
  
  # Process each chunk file sequentially
  results <- lapply(chunk_files, process_chunk)
  
  # Combine all annual rasters into a single stack for the selected metrics
  combine_rasters <- function(raster_list) {
    raster_stack <- rast()
    for (year in names(raster_list)) {
      if (!is.null(raster_list[[year]])) {
        raster_stack <- c(raster_stack, raster_list[[year]])
      }
    }
    return(raster_stack)
  }
  
  persistence_stack <- NULL
  duration_stack <- NULL
  concentration_stack <- NULL
  
  if ("persistence" %in% report_metrics) {
    persistence_stack <- combine_rasters(do.call(c, lapply(results, `[[`, "persistence")))
    if (use_mask) persistence_stack <- mask(persistence_stack, study_area)
  }
  
  if ("duration" %in% report_metrics) {
    duration_stack <- combine_rasters(do.call(c, lapply(results, `[[`, "duration")))
    if (use_mask) duration_stack <- mask(duration_stack, study_area)
  }
  
  if ("concentration" %in% report_metrics) {
    concentration_stack <- combine_rasters(do.call(c, lapply(results, `[[`, "concentration")))
    if (use_mask) concentration_stack <- mask(concentration_stack, study_area)
  }
  
  # Function to apply linear model to each pixel and return slope and p-value
  apply_lm <- function(y) {
    if (all(is.na(y))) {
      return(c(NA, NA))
    } else {
      years <- start_year:end_year
      fit <- tryCatch(lm(y ~ years), error = function(e) return(NULL))
      if (is.null(fit)) {
        return(c(NA, NA))
      } else {
        slope <- coef(fit)[2]  # Slope of the linear model
        p_value <- summary(fit)$coefficients[2, 4]  # p-value for the slope
        return(c(slope, p_value))
      }
    }
  }
  
  persistence_trend <- NULL
  duration_trend <- NULL
  concentration_trend <- NULL
  
  if (!is.null(persistence_stack)) persistence_trend <- app(persistence_stack, apply_lm)
  if (!is.null(duration_stack)) duration_trend <- app(duration_stack, apply_lm)
  if (!is.null(concentration_stack)) concentration_trend <- app(concentration_stack, apply_lm)
  
  # Extract slopes and p-values from the results
  persistence_slope <- NULL
  persistence_pvalue <- NULL
  duration_slope <- NULL
  duration_pvalue <- NULL
  concentration_slope <- NULL
  concentration_pvalue <- NULL
  
  if (!is.null(persistence_trend)) {
    persistence_slope <- persistence_trend[[1]]
    persistence_pvalue <- persistence_trend[[2]]
  }
  
  if (!is.null(duration_trend)) {
    duration_slope <- duration_trend[[1]]
    duration_pvalue <- duration_trend[[2]]
  }
  
  if (!is.null(concentration_trend)) {
    concentration_slope <- concentration_trend[[1]]
    concentration_pvalue <- concentration_trend[[2]]
  }
  
  # Verify the integrity of the data
  if (!is.null(persistence_slope)) {
    print("Summary of persistence slope layer:")
    print(summary(persistence_slope))
    print("Summary of persistence p-value layer:")
    print(summary(persistence_pvalue))
  }
  
  if (!is.null(duration_slope)) {
    print("Summary of duration slope layer:")
    print(summary(duration_slope))
    print("Summary of duration p-value layer:")
    print(summary(duration_pvalue))
  }
  
  if (!is.null(concentration_slope)) {
    print("Summary of concentration slope layer:")
    print(summary(concentration_slope))
    print("Summary of concentration p-value layer:")
    print(summary(concentration_pvalue))
  }
  
  # Define a significance level
  alpha <- 0.05
  
  # Calculate statistics for persistence
  if (!is.null(persistence_slope)) {
    persistence_values <- values(persistence_slope)
    persistence_pvalues <- values(persistence_pvalue)
    persistence_data <- cbind(persistence_values, persistence_pvalues)
    persistence_data <- persistence_data[complete.cases(persistence_data), ]
    persistence_values <- persistence_data[, 1]
    persistence_pvalues <- persistence_data[, 2]
    persistence_mean <- mean(persistence_values)
    persistence_median <- median(persistence_values)
    persistence_min <- min(persistence_values)
    persistence_max <- max(persistence_values)
    persistence_significant <- sum(persistence_pvalues <= alpha) / length(persistence_pvalues) * 100
    persistence_positive <- sum(persistence_pvalues <= alpha & persistence_values > 0) / length(persistence_pvalues) * 100
    persistence_negative <- sum(persistence_pvalues <= alpha & persistence_values < 0) / length(persistence_pvalues) * 100
    persistence_nonsignificant <- 100 - persistence_significant
  }
  
  # Calculate statistics for duration
  if (!is.null(duration_slope)) {
    duration_values <- values(duration_slope)
    duration_pvalues <- values(duration_pvalue)
    duration_data <- cbind(duration_values, duration_pvalues)
    duration_data <- duration_data[complete.cases(duration_data), ]
    duration_values <- duration_data[, 1]
    duration_pvalues <- duration_data[, 2]
    duration_mean <- mean(duration_values)
    duration_median <- median(duration_values)
    duration_min <- min(duration_values)
    duration_max <- max(duration_values)
    duration_significant <- sum(duration_pvalues <= alpha) / length(duration_pvalues) * 100
    duration_positive <- sum(duration_pvalues <= alpha & duration_values > 0) / length(duration_pvalues) * 100
    duration_negative <- sum(duration_pvalues <= alpha & duration_values < 0) / length(duration_pvalues) * 100
    duration_nonsignificant <- 100 - duration_significant
  }
  
  # Calculate statistics for concentration
  if (!is.null(concentration_slope)) {
    concentration_values <- values(concentration_slope)
    concentration_pvalues <- values(concentration_pvalue)
    concentration_data <- cbind(concentration_values, concentration_pvalues)
    concentration_data <- concentration_data[complete.cases(concentration_data), ]
    concentration_values <- concentration_data[, 1]
    concentration_pvalues <- concentration_data[, 2]
    concentration_mean <- mean(concentration_values)
    concentration_median <- median(concentration_values)
    concentration_min <- min(concentration_values)
    concentration_max <- max(concentration_values)
    concentration_significant <- sum(concentration_pvalues <= alpha) / length(concentration_pvalues) * 100
    concentration_positive <- sum(concentration_pvalues <= alpha & concentration_values > 0) / length(concentration_pvalues) * 100
    concentration_negative <- sum(concentration_pvalues <= alpha & concentration_values < 0) / length(concentration_pvalues) * 100
    concentration_nonsignificant <- 100 - concentration_significant
  }
  
  # Print summary statistics for reporting
  if (!is.null(persistence_slope)) {
    persistence_summary <- data.frame(
      Metric = c("Mean slope", "Median slope", "Minimum slope", "Maximum slope", "% Significant", "% Positive", "% Negative", "% Non-significant"),
      Persistence = c(persistence_mean, persistence_median, persistence_min, persistence_max, persistence_significant, persistence_positive, persistence_negative, persistence_nonsignificant)
    )
    print(persistence_summary)
  }
  
  if (!is.null(duration_slope)) {
    duration_summary <- data.frame(
      Metric = c("Mean slope", "Median slope", "Minimum slope", "Maximum slope", "% Significant", "% Positive", "% Negative", "% Non-significant"),
      Duration = c(duration_mean, duration_median, duration_min, duration_max, duration_significant, duration_positive, duration_negative, duration_nonsignificant)
    )
    print(duration_summary)
  }
  
  if (!is.null(concentration_slope)) {
    concentration_summary <- data.frame(
      Metric = c("Mean slope", "Median slope", "Minimum slope", "Maximum slope", "% Significant", "% Positive", "% Negative", "% Non-significant"),
      Concentration = c(concentration_mean, concentration_median, concentration_min, concentration_max, concentration_significant, concentration_positive, concentration_negative, concentration_nonsignificant)
    )
    print(concentration_summary)
  }
  
  # Define a color palette for trends
  trend_colors <- colorRampPalette(c("red", "white", "blue"))
  
  # Set up multi-panel plotting environment
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  
  if (!is.null(persistence_slope)) {
    plot(persistence_slope, col=trend_colors(100), main="Trend of Persistence Over Time", legend=TRUE)
  }
  
  if (!is.null(duration_slope)) {
    plot(duration_slope, col=trend_colors(100), main="Trend of Duration Over Time", legend=TRUE)
  }
  
  if (!is.null(concentration_slope)) {
    plot(concentration_slope, col=trend_colors(100), main="Trend of Concentration Over Time", legend=TRUE)
  }
  
  # Classify the trends based on significance and directionality using ifel function
  if (!is.null(persistence_pvalue)) {
    persistence_significance_direction <- ifel(persistence_pvalue <= alpha & persistence_slope > 0, 1,
                                               ifel(persistence_pvalue <= alpha & persistence_slope < 0, -1, 0))
    plot(persistence_significance_direction, col=c("red", "white", "blue"), main="Significance and Directionality of Persistence Trends", legend=TRUE)
  }
  
  if (!is.null(duration_pvalue)) {
    duration_significance_direction <- ifel(duration_pvalue <= alpha & duration_slope > 0, 1,
                                            ifel(duration_pvalue <= alpha & duration_slope < 0, -1, 0))
    plot(duration_significance_direction, col=c("red", "white", "blue"), main="Significance and Directionality of Duration Trends", legend=TRUE)
  }
  
  if (!is.null(concentration_pvalue)) {
    concentration_significance_direction <- ifel(concentration_pvalue <= alpha & concentration_slope > 0, 1,
                                                 ifel(concentration_pvalue <= alpha & concentration_slope < 0, -1, 0))
    plot(concentration_significance_direction, col=c("red", "white", "blue"), main="Significance and Directionality of Concentration Trends", legend=TRUE)
  }
}

# Example usage
analyze_sea_ice(start_year = 2014, end_year = 2023, use_mask = FALSE, report_metrics = c("concentration", "duration", "persistence"))
