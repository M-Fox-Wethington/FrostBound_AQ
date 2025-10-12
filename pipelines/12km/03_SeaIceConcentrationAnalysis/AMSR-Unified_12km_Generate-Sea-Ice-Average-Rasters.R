library(terra)
library(sf)

# Load the raster stack
r_fp <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/AMSR-Unified_SIC_Full-Catalog_Time_Series_Subregions.nc"
r_stack <- rast(r_fp)

# Apply SIC > 15 filter directly on the stack
r_filtered <- ifel(r_stack > 15, r_stack, NA)

# Ensure the time dimension is available and correctly formatted
if (is.na(time(r_filtered)[1])) {
  message("Time dimension not properly set. Setting time dimension...")
  start_date <- as.Date("2012-07-03")  # Adjust according to your actual data
  time_values <- seq(start_date, length.out = nlyr(r_stack), by = "day")
  time(r_filtered) <- time_values
}

# Define winter months for the Southern Hemisphere
winter_months <- c("06", "07", "08", "09")

# Prepare containers for results
winter_monthly_means <- list()
annual_winter_means <- list()

# Extract dates, years, and months
dates <- as.Date(time(r_filtered))
years <- format(dates, "%Y")
months <- format(dates, "%m")

# Aggregate by winter months for each year
for (year in unique(years)) {
  # Calculate mean for each winter month individually
  for (month in winter_months) {
    idx <- which(years == year & months == month)
    if (length(idx) > 0) {
      winter_monthly_means[[paste(year, month, sep="-")]] <- mean(r_filtered[[idx]], na.rm = TRUE)
    }
  }
  
  # Calculate mean for the entire winter season each year
  winter_idx <- which(years == year & months %in% winter_months)
  if (length(winter_idx) > 0) {
    annual_winter_means[[year]] <- mean(r_filtered[[winter_idx]], na.rm = TRUE)
  }
}

# Convert lists of rasters to SpatRasters
winter_monthly_mean_raster <- rast(winter_monthly_means)
annual_winter_mean_raster <- rast(annual_winter_means)

# Name layers by year-month for monthly and by year for seasonal
names(winter_monthly_mean_raster) <- names(winter_monthly_means)
names(annual_winter_mean_raster) <- names(annual_winter_means)


# Paths for output NetCDF files
month_raststack_output_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/substack-extents/AMSR-Unified_winter_monthly_mean_extent.nc"
annual_raststack_output_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/12km_AMSR-Unified/stack/substack-extents/AMSR-Unified_annual_winter_mean_extent.nc"

# Write results to NetCDF files using writeCDF
writeCDF(winter_monthly_mean_raster, month_raststack_output_path, varname="Winter_Monthly_Mean_Ice_Extent", overwrite=TRUE)
writeCDF(annual_winter_mean_raster, annual_raststack_output_path, varname="Winter_Annual_Mean_Ice_Extent", overwrite=TRUE)
