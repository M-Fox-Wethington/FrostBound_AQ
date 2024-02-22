# Load required libraries
# terra: For spatial data analysis and raster manipulation
# rhdf5: For reading HDF5 files, a format commonly used for storing complex scientific data
# sf: For handling spatial vector data (e.g., points, lines, polygons)
library(terra)
library(rhdf5)
library(sf)

# Set the file path for the HDF5 file containing Sea Ice Concentration data
h5_file <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/178568954/AMSR_U2_L3_SeaIce12km_B04_20120704.he5"

# Read Sea Ice Concentration data from the specified HDF5 dataset path
SeaIce <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/Data Fields/SI_12km_SH_ICECON_DAY")

# Read and display detailed metadata about the dataset
# This metadata can include information about the data collection, processing, and file structure
Core_Metadata <- h5read(h5_file, "/HDFEOS INFORMATION/CoreMetadata.0")
print(Core_Metadata)

Struct_Metadata <- h5read(h5_file, "/HDFEOS INFORMATION/StructMetadata.0")
print(Struct_Metadata)
str(Struct_Metadata)  # Display the structure of the Struct_Metadata for understanding its hierarchy

# Read spatial dimension data which represents the size of the grid in the X and Y directions
XDim <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/XDim")
print(XDim)
YDim <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/YDim")
print(YDim)

# Read latitude and longitude data for georeferencing the Sea Ice data
# These datasets provide the geographic coordinates for each cell in the Sea Ice dataset
lat <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/lat")
lon <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/lon")

# Get the dimensions of the latitude (and longitude) data to understand the grid size
dims <- dim(lat)

# Flatten the Sea Ice, latitude, and longitude matrices into vectors to facilitate creating point data
SeaIce_vec <- as.vector(SeaIce)
lat_vec <- as.vector(lat)
lon_vec <- as.vector(lon)

# Combine the flattened vectors into a data frame and then convert it into an 'sf' (Simple Features) object
# This creates a spatially-aware dataset that associates each Sea Ice Concentration value with a geographic location
sic_points <- data.frame(lon = lon_vec, lat = lat_vec, SIC = SeaIce_vec)
sic_points_sf <- st_as_sf(sic_points, coords = c("lon", "lat"), crs = 4326, agr = "constant")

# Transform the geographic coordinates from WGS84 (EPSG:4326) to the Antarctic Polar Stereographic South projection (EPSG:3412)
# This projection is appropriate for spatial analysis related to the Antarctic region
sic_points_sf_transformed <- st_transform(sic_points_sf, crs = "EPSG:3412")

# Create a raster template with dimensions and extent matching the transformed spatial points
# This raster will be used to rasterize the point data, assigning each point's SIC value to the appropriate raster cell
r_corrected <- rast(nrows=632, ncols=664, xmin=-3943750, xmax=3943750, ymin=-3943750, ymax=4343750, crs="EPSG:3412")

# Rasterize the transformed point data onto the raster template
# The 'mean' function is used to aggregate multiple points that may fall into a single raster cell
r_sic <- rasterize(sic_points_sf_transformed, r_corrected, field="SIC", fun=mean)

# Visualize the resulting rasterized Sea Ice Concentration data
plot(r_sic)

# Save the raster data to a GeoTIFF file, overwriting the file if it already exists
# GeoTIFF is a widely used format for storing raster data with embedded georeferencing information
writeRaster(r_sic, filename="D:/Manuscripts_localData/FrostBound_AQ/Datasets/Georeferenced_SeaIce_Concentration.tif", overwrite = TRUE)
