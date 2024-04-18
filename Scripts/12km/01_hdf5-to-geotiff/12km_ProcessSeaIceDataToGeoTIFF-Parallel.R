# Load required libraries
library(terra)   # For spatial data analysis and raster manipulation
library(rhdf5)   # For reading HDF5 files
library(sf)      # For handling spatial vector data
library(doParallel)
library(foreach)

# Ensure output directory exists
output_dir <- "D:/D-Downloads/AMSRE-Unified_Catalog/processed"

# output_dir <-  "D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-Data/Processed"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Define a function to process a single HDF5 file
process_hdf5_file <- function(h5_file, input_dir, output_dir) {
  # Extract the date from the file path using regex, assuming YYYYMMDD format before ".he5"
  date_pattern <- ".*(\\d{4})(\\d{2})(\\d{2})\\.he5$"
  date <- sub(date_pattern, "\\1-\\2-\\3", basename(h5_file))
  
  # Read Sea Ice Concentration data and apply scale factor
  SeaIce_scaled <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/Data Fields/SI_12km_SH_ICECON_DAY")
  
  
  # Set value ranges based on documentation
  SeaIce_scaled[SeaIce_scaled == 110] <- NA  # Missing data
  SeaIce_scaled[SeaIce_scaled == 120] <- NA  # Land
  
  # Read spatial dimension data for georeferencing
  lat <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/lat")
  lon <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/lon")
  
  # Flatten the matrices to vectors
  SeaIce_vec <- as.vector(SeaIce_scaled)
  lat_vec <- as.vector(lat)
  lon_vec <- as.vector(lon)
  
  # Combine vectors into a data frame and convert to an 'sf' object
  sic_points <- data.frame(lon = lon_vec, lat = lat_vec, SIC = SeaIce_vec)
  sic_points_sf <- st_as_sf(sic_points, coords = c("lon", "lat"), crs = 4326, agr = "constant")
  
  # Transform geographic coordinates to the Antarctic Polar Stereographic South projection (EPSG:3412)
  sic_points_sf_transformed <- st_transform(sic_points_sf, crs = "EPSG:3412")
  
  # Create a raster template with dimensions and extent matching the spatial points
  r_corrected <- rast(nrows=664, ncols=632, xmin=-3950000, xmax=3950000, ymin=-3950000, ymax=4350000, crs="EPSG:3412")
  
  # Rasterize the transformed point data onto the raster template
  r_sic <- rasterize(sic_points_sf_transformed, r_corrected, field="SIC", fun=mean)
  
  # Define the output filename using the output directory and date
  output_filename <- paste0(output_dir, "/NSIDC-AU_SI12_SeaIce_12km_", gsub("-", "", date), ".tif")
  
  # Save the GeoTIFF in the output directory
  writeRaster(r_sic, filename=output_filename, overwrite = TRUE)
  
  print(paste0("Raster:", output_filename, " exported"))
  
  # Clear memory of large objects no longer needed
  rm(date, date_pattern, r_sic, sic_points, sic_points_sf, sic_points_sf_transformed, SeaIce_scaled, SeaIce_vec, lat, lat_vec, lon, lon_vec)
  
  # Explicitly call garbage collection to free up memory space
  gc()
}


# Setup parallel environment
numCores <- detectCores()
cl <- makeCluster(numCores - 1)
registerDoParallel(cl)# Load required libraries
library(terra)   # For spatial data analysis and raster manipulation
library(rhdf5)   # For reading HDF5 files
library(sf)      # For handling spatial vector data
library(doParallel)
library(foreach)

# Ensure output directory exists
output_dir <- "D:/D-Downloads/AMSR-Unified_Catalog/processed"

# output_dir <-  "D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-Data/Processed"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Define a function to process a single HDF5 file
process_hdf5_file <- function(h5_file, input_dir, output_dir) {
  # Extract the date from the file path using regex, assuming YYYYMMDD format before ".he5"
  date_pattern <- ".*(\\d{4})(\\d{2})(\\d{2})\\.he5$"
  date <- sub(date_pattern, "\\1-\\2-\\3", basename(h5_file))
  
  # Read Sea Ice Concentration data and apply scale factor
  SeaIce_scaled <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/Data Fields/SI_12km_SH_ICECON_DAY")
  
  
  # Set value ranges based on documentation
  SeaIce_scaled[SeaIce_scaled == 110] <- NA  # Missing data
  SeaIce_scaled[SeaIce_scaled == 120] <- NA  # Land
  
  # Read spatial dimension data for georeferencing
  lat <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/lat")
  lon <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/lon")
  
  # Flatten the matrices to vectors
  SeaIce_vec <- as.vector(SeaIce_scaled)
  lat_vec <- as.vector(lat)
  lon_vec <- as.vector(lon)
  
  # Combine vectors into a data frame and convert to an 'sf' object
  sic_points <- data.frame(lon = lon_vec, lat = lat_vec, SIC = SeaIce_vec)
  sic_points_sf <- st_as_sf(sic_points, coords = c("lon", "lat"), crs = 4326, agr = "constant")
  
  # Transform geographic coordinates to the Antarctic Polar Stereographic South projection (EPSG:3412)
  sic_points_sf_transformed <- st_transform(sic_points_sf, crs = "EPSG:3412")
  
  # Create a raster template with dimensions and extent matching the spatial points
  r_corrected <- rast(nrows=664, ncols=632, xmin=-3950000, xmax=3950000, ymin=-3950000, ymax=4350000, crs="EPSG:3412")
  
  # Rasterize the transformed point data onto the raster template
  r_sic <- rasterize(sic_points_sf_transformed, r_corrected, field="SIC", fun=mean)
  
  # Define the output filename using the output directory and date
  output_filename <- paste0(output_dir, "/NSIDC-AU_SI12_SeaIce_12km_", gsub("-", "", date), ".tif")
  
  # Save the GeoTIFF in the output directory
  writeRaster(r_sic, filename=output_filename, overwrite = TRUE)
  
  print(paste0("Raster:", output_filename, " exported"))
  
  # Clear memory of large objects no longer needed
  rm(date, date_pattern, r_sic, sic_points, sic_points_sf, sic_points_sf_transformed, SeaIce_scaled, SeaIce_vec, lat, lat_vec, lon, lon_vec)
  
  # Explicitly call garbage collection to free up memory space
  gc()
}


# Setup parallel environment
numCores <- detectCores()
cl <- makeCluster(numCores - 1)
registerDoParallel(cl)

# Export necessary objects to the cluster
clusterExport(cl, varlist = c("process_hdf5_file", "input_dir", "output_dir"))

# Directory containing HDF5 files
input_dir <- "D:/D-Downloads/AMSR-Unified_Catalog/staged"

# List all HDF5 files in the directory
hdf5_files <- list.files(input_dir, pattern = "\\.he5$", full.names = TRUE)

# Apply the processing function to each HDF5 file in parallel
results <- foreach(h5_file = hdf5_files, .packages = c("terra", "rhdf5", "sf")) %dopar% {
  process_hdf5_file(h5_file, input_dir, output_dir)
}

# Stop the cluster
stopCluster(cl)

# Export necessary objects to the cluster
clusterExport(cl, varlist = c("process_hdf5_file", "input_dir", "output_dir"))

# Directory containing HDF5 files
input_dir <- "D:/D-Downloads/AMSR-Unified_Catalog/staged"

# List all HDF5 files in the directory
hdf5_files <- list.files(input_dir, pattern = "\\.he5$", full.names = TRUE)

# Apply the processing function to each HDF5 file in parallel
results <- foreach(h5_file = hdf5_files, .packages = c("terra", "rhdf5", "sf")) %dopar% {
  process_hdf5_file(h5_file, input_dir, output_dir)
}

# Stop the cluster
stopCluster(cl)