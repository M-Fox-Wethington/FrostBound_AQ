library(terra)
library(lubridate)

# Load the NSIDC sea ice concentration data
nsidc <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/stack/NSIDC_25km_Study-Area.nc")

# Define a function to count days with ice cover > 15% for each year
count_ice_days <- function(x) {
  ice_days <- sum(x > 15, na.rm = TRUE)
  return(ice_days)
}

# Extract the years from the time attribute using lubridate
years <- unique(year(time(nsidc)))

# Initialize an empty list to store the annual ice day count rasters
annual_ice_days <- list()

# Loop over each year, subset the raster stack, and apply the counting function
for (year in years) {
  year_rasters <- subset(nsidc, year(time(nsidc)) == year)
  # Apply count_ice_days function across layers (days) at each pixel
  ice_days_raster <- app(year_rasters, count_ice_days)  # No need for 'along=3', app inherently operates cell-wise
  names(ice_days_raster) <- paste("Ice_Days", year, sep="_")
  annual_ice_days[[year]] <- ice_days_raster
  print(paste("Processed year:", year))
}

# Combine the annual rasters into a single stack
annual_ice_days_stack <- rast(annual_ice_days)

# Save the resulting stack to disk
# writeRaster(annual_ice_days_stack, "D:/Manuscripts_localData/FrostBound_AQ/Output/Annual_Ice_Days_Stack.tif", overwrite=TRUE)

# Optionally, plot one of the layers to visualize
plot(annual_ice_days_stack$Ice_Days_2020, main="Days with >15% Ice Cover in 2020")


plot(annual_ice_days_stack$Ice_Days_1980)
