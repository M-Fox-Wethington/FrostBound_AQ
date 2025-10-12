# Load necessary libraries
library(terra)
library(dplyr)
library(lubridate)
library(nlme)
library(ggplot2)
library(grDevices)  # For colorRampPalette
library(gridExtra)  # For arranging plots side by side

# Helper function to calculate duration
calculate_duration <- function(x, threshold) {
  if (all(is.na(x))) return(NA)
  rle_result <- rle(x > threshold)
  max_duration <- ifelse(any(rle_result$values), max(rle_result$lengths[rle_result$values]), 0)
  return(max_duration)
}

# Helper function to calculate linear model for each pixel
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

# Helper function to calculate statistics
calculate_statistics <- function(values, pvalues, alpha) {
  data <- cbind(values, pvalues)
  data <- data[complete.cases(data), ]
  values <- data[, 1]
  pvalues <- data[, 2]
  
  mean_value <- mean(values)
  median_value <- median(values)
  min_value <- min(values)
  max_value <- max(values)
  significant <- sum(pvalues <= alpha) / length(pvalues) * 100
  positive <- sum(pvalues <= alpha & values > 0) / length(pvalues) * 100
  negative <- sum(pvalues <= alpha & values < 0) / length(pvalues) * 100
  nonsignificant <- 100 - significant
  
  return(list(mean = mean_value, median = median_value, min = min_value, max = max_value, 
              significant = significant, positive = positive, negative = negative, nonsignificant = nonsignificant))
}

# Helper function to print summary statistics
print_summary <- function(metric_name, statistics) {
  summary_df <- data.frame(
    Metric = c("Mean slope", "Median slope", "Minimum slope", "Maximum slope", "% Significant", "% Positive", "% Negative", "% Non-significant"),
    Value = c(statistics$mean, statistics$median, statistics$min, statistics$max, statistics$significant, statistics$positive, statistics$negative, statistics$nonsignificant)
  )
  colnames(summary_df)[2] <- metric_name
  print(summary_df)
}

# Helper function to classify and plot trends
plot_trends <- function(slope, pvalue, alpha, title) {
  significant_cells <- pvalue <= alpha
  significant_mask <- ifel(significant_cells, 1, NA)
  positive_trend <- ifel(significant_mask == 1 & slope > 0, 1, NA)
  negative_trend <- ifel(significant_mask == 1 & slope < 0, 2, NA)
  na_mask <- ifel(is.na(slope), 3, NA)
  combined_trend <- positive_trend
  combined_trend[!is.na(negative_trend)] <- 2
  combined_trend[!is.na(na_mask)] <- 3
  color_table <- data.frame(value = c(1, 2, 3), col = c("blue", "red", "black"))
  coltab(combined_trend) <- color_table
  plot(combined_trend, main = title)
  legend("topright", legend = c("Positive Trend", "Negative Trend", "NA Cells"), fill = color_table$col)
}

# Define the main function
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
        persistence <- app(nsidc_year, function(x) sum(ifelse(x < ice_threshold, 0, 1), na.rm = FALSE))
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
  
  # Calculate and print statistics for each metric
  if (!is.null(persistence_slope)) {
    persistence_statistics <- calculate_statistics(values(persistence_slope), values(persistence_pvalue), alpha)
    print_summary("Persistence", persistence_statistics)
  }
  
  if (!is.null(duration_slope)) {
    duration_statistics <- calculate_statistics(values(duration_slope), values(duration_pvalue), alpha)
    print_summary("Duration", duration_statistics)
  }
  
  if (!is.null(concentration_slope)) {
    concentration_statistics <- calculate_statistics(values(concentration_slope), values(concentration_pvalue), alpha)
    print_summary("Concentration", concentration_statistics)
  }
  
  # Define a color palette for trends
  trend_colors <- colorRampPalette(c("red", "white", "blue"))
  
  # Set up multi-panel plotting environment
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  
  if (!is.null(persistence_slope)) {
    plot(persistence_slope, col=trend_colors(100), main="Trend of Persistence Over Time", legend=TRUE, colNA="black")
  }
  
  if (!is.null(duration_slope)) {
    plot(duration_slope, col=trend_colors(100), main="Trend of Duration Over Time", legend=TRUE, colNA="black")
  }
  
  if (!is.null(concentration_slope)) {
    plot(concentration_slope, col=trend_colors(100), main="Trend of Concentration Over Time", legend=TRUE, colNA="black")
  }
  
  # Classify and plot trends based on significance and directionality
  if (!is.null(persistence_pvalue)) {
    plot_trends(persistence_slope, persistence_pvalue, alpha, "Significance and Directionality of Persistence Trends")
  }
  
  if (!is.null(duration_pvalue)) {
    plot_trends(duration_slope, duration_pvalue, alpha, "Significance and Directionality of Duration Trends")
  }
  
  if (!is.null(concentration_pvalue)) {
    plot_trends(concentration_slope, concentration_pvalue, alpha, "Significance and Directionality of Concentration Trends")
  }
}

# Example usage
analyze_sea_ice(start_year = 2014, end_year = 2023, use_mask = FALSE, report_metrics = c("concentration", "duration", "persistence"))
