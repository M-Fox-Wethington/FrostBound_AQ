library(terra)

# Define a function to extract the date and return as Date object
extractDate <- function(filename) {
  base_name <- basename(filename)
  
  # capture the date directly following '_12km_'
  date_string <- sub(".*_12km_(\\d{8}).*", "\\1", base_name)
  
  # Attempt to convert to Date
  tryCatch({
    as.Date(date_string, format="%Y%m%d")
  }, warning = function(w) {
    warning("Date conversion failed for: ", base_name)
    NA
  }, error = function(e) {
    warning("Date conversion failed for: ", base_name)
    NA
  })
}



# Load raster files
raster_dir <- "D:/D-Downloads/AMSR-Unified_Catalog/processed"
raster_files <- list.files(path = raster_dir, pattern = "\\.tif$", full.names = TRUE)
raster_files <- raster_files[order(sapply(raster_files, extractDate))]

# Create a raster stack
raster_stack <- rast(lapply(raster_files, rast))
names(raster_stack) <- sapply(raster_files, function(x) format(extractDate(x), "%Y-%m-%d"))

# Assigning time values
dates <- as.Date(sapply(raster_files, extractDate))
time(raster_stack) <- dates

# Ensure CRS matches, clip and mask
region_shp_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/CCAMLR_MPA_Planning_Domains_WAP.shp"

region <- vect(region_shp_path)
sample_raster <- raster_stack[[1]]
if (crs(region) != crs(sample_raster)) {
  region <- project(region, crs(sample_raster))
}
clipped_stack <- crop(raster_stack, ext(region))
clipped_stack <- mask(clipped_stack, region)

# Save and reload the raster stack
output_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-Unified_12km/AMSR-Unified_SIC_Time_Series.nc"
writeCDF(clipped_stack, filename = output_path, overwrite = TRUE)

loaded_stack <- rast(output_path)
print(loaded_stack)
print(time(loaded_stack))


# Assuming 'loaded_stack' is your loaded SpatRaster with supposed time data
time_data <- time(loaded_stack)
if (length(time_data) > 0 && !all(is.na(time_data))) {
  print("Time data are present and appear correctly set.")
  print(head(time_data, 20))  # Print the first 20 time entries to verify
} else {
  print("No valid time data found or all entries are NA.")
}
