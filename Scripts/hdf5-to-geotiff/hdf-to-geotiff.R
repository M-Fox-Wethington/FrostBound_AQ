# Load required libraries
library(terra)   # For spatial data analysis and raster manipulation
library(rhdf5)   # For reading HDF5 files, a format commonly used for storing complex scientific data
library(sf)      # For handling spatial vector data (e.g., points, lines, polygons)

# Set the file path for the HDF5 file containing Sea Ice Concentration data
h5_file <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/178568954/AMSR_U2_L3_SeaIce12km_B04_20120704.he5"

# Extract the date from the file path using regex, assuming YYYYMMDD format before ".he5"
date_pattern <- ".*(\\d{4})(\\d{2})(\\d{2})\\.he5$"
date <- sub(date_pattern, "\\1-\\2-\\3", basename(h5_file))

# Split the date to extract the year and month parts separately
date_parts <- strsplit(date, "-")[[1]]
year <- date_parts[1]
month <- date_parts[2]

# Create a new subdirectory path based on the year and month
new_subdir <- paste0("D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-Data/Processed/", year, "/", month, "/")

# Ensure the new subdirectory exists, if not, create it
if (!dir.exists(new_subdir)) {
  dir.create(new_subdir, recursive = TRUE)
}

# Read Sea Ice Concentration data from the HDF5 dataset
SeaIce <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/Data Fields/SI_12km_SH_ICECON_DAY")

# Optionally, read and display detailed metadata about the dataset
Core_Metadata <- h5read(h5_file, "/HDFEOS INFORMATION/CoreMetadata.0")
Struct_Metadata <- h5read(h5_file, "/HDFEOS INFORMATION/StructMetadata.0")

# Read spatial dimension data (XDim and YDim) for grid size
XDim <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/XDim")
YDim <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/YDim")

# Read latitude and longitude data for georeferencing
lat <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/lat")
lon <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/lon")

# Flatten the Sea Ice, latitude, and longitude matrices into vectors
SeaIce_vec <- as.vector(SeaIce)
lat_vec <- as.vector(lat)
lon_vec <- as.vector(lon)

# Combine vectors into a data frame and convert to an 'sf' object
sic_points <- data.frame(lon = lon_vec, lat = lat_vec, SIC = SeaIce_vec)
sic_points_sf <- st_as_sf(sic_points, coords = c("lon", "lat"), crs = 4326, agr = "constant")

# Transform geographic coordinates to the Antarctic Polar Stereographic South projection (EPSG:3412)
sic_points_sf_transformed <- st_transform(sic_points_sf, crs = "EPSG:3412")

# Create a raster template with dimensions and extent matching the spatial points
r_corrected <- rast(nrows=632, ncols=664, xmin=-3943750, xmax=3943750, ymin=-3943750, ymax=4343750, crs="EPSG:3412")

# Rasterize the transformed point data onto the raster template
r_sic <- rasterize(sic_points_sf_transformed, r_corrected, field="SIC", fun=mean)

# Visualize the rasterized Sea Ice Concentration data
plot(r_sic)

# Define the output filename using the new subdirectory and date
output_filename <- paste0(new_subdir, "SeaIce_Concentration_", gsub("-", "", date), ".tif")

# Save the GeoTIFF in the newly created subdirectory
writeRaster(r_sic, filename=output_filename, overwrite = TRUE)
