# Load required libraries
library(terra)   # For spatial data analysis and raster manipulation
library(sf)      # For handling spatial vector data (e.g., points, lines, polygons)

# Define a function to process a single NetCDF file
process_netcdf_file <- function(nc_file, input_dir, output_dir) {
  # Extract the date from the file path using regex, assuming YYYYMMDD format
  date_pattern <- ".*_(\\d{4})(\\d{2})(\\d{2})_v\\d+\\.\\d+\\.nc$"
  date <- sub(date_pattern, "\\1-\\2-\\3", basename(nc_file))
  
  # Read the sea ice concentration data from the NetCDF file
  sic_raster <- rast(nc_file, subds="Sea_ice_concentration") # Adjust "Sea_ice_concentration" based on the actual variable name
  
  # Optionally, transform the projection if it's not already in the desired CRS
  # sic_raster_transformed <- project(sic_raster, "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  # Note: Uncomment and adjust the above line if a projection transformation is needed.
  
  # Define the output filename using the specified output directory and date
  output_filename <- sprintf("%sNSIDC0079_SEAICE_PS_S25km_%s.tif", output_dir, gsub("-", "", date))
  
  # Save the GeoTIFF in the designated output directory
  writeRaster(sic_raster, filename=output_filename, overwrite=TRUE)
  
  print(paste0("Raster:", output_filename, " exported"))
  
  # Clear memory of large objects no longer needed
  rm(date, output_filename, sic_raster)
  
  # Explicitly call garbage collection to free up memory space
  gc()
}

# Directory containing NetCDF files
input_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/SSMI/tmp"
# Single output directory for all TIFF files
output_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/SSMI/tmp/output/"

# Ensure the output directory exists, if not, create it
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# List all NetCDF files in the directory
netcdf_files <- list.files(input_dir, pattern = "\\.nc$", full.names = TRUE)

# Apply the processing function to each NetCDF file
lapply(netcdf_files, function(file) process_netcdf_file(file, input_dir, output_dir))
