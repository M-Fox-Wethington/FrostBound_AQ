library(terra)
library(lubridate)

process_winter_sic <- function(raster_path, threshold = 15) {
  # Load the raster stack
  r_stack <- rast(raster_path)
  
  # Define winter months for the Southern Hemisphere
  winter_months <- c("06", "07", "08", "09")
  
  # Check if the raster has time data and process accordingly
  if (!is.null(time(r_stack))) {
    months <- format(time(r_stack), "%m")
    winter_idx <- which(months %in% winter_months)
    r_winter <- r_stack[[winter_idx]]
  } else {
    stop("Raster does not have time data")
  }
  
  # Apply a threshold to remove low values (less meaningful for ice concentration)
  r_winter[r_winter <= threshold] <- NA
  
  # Prepare containers for results
  winter_monthly_means <- list()
  annual_winter_means <- list()
  
  # Re-extract formatted years and months for subsetting
  years <- format(time(r_winter), "%Y")
  months <- format(time(r_winter), "%m")
  
  # Aggregate by winter months for each year
  for (year in unique(years)) {
    for (month in winter_months) {
      idx <- which(years == year & months == month)
      if (length(idx) > 0) {
        winter_monthly_means[[paste(year, month, sep="-")]] <- mean(r_winter[[idx]], na.rm = TRUE)
      }
    }
    winter_idx <- which(years == year & months %in% winter_months)
    if (length(winter_idx) > 0) {
      annual_winter_means[[year]] <- mean(r_winter[[winter_idx]], na.rm = TRUE)
    }
  }
  
  # Convert lists of rasters to SpatRasters
  winter_monthly_mean_raster <- rast(winter_monthly_means)
  annual_winter_mean_raster <- rast(annual_winter_means)
  
  # Optionally, plot or save results
  plot(winter_monthly_mean_raster, main=paste("Monthly Winter Mean SIC for", basename(raster_path)))
  plot(annual_winter_mean_raster, main=paste("Annual Winter Mean SIC for", basename(raster_path)))
  
  # Return the processed data
  return(list(monthly=winter_monthly_mean_raster, annual=annual_winter_mean_raster))
}

# Use the function on AMSR and NSIDC datasets
amsr_results <- process_winter_sic("D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-Unified_12km/stack/AMSR-Unified_12km_Harmonization.nc")
nsidc_results <- process_winter_sic("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/25km_Sea-Ice-Index/stack/NSIDC-Sea-Ice-Index_25km_Winter_Harmonization.nc")
