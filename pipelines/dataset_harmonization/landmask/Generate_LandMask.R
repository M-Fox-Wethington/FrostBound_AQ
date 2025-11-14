# Load necessary libraries
library(terra)

# Load the NetCDF file
nc_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/Datasets/25km_Sea-Ice-Index/processed_clipped/S_20220307_concentration_v3.0.tif"
nsidc_raster <- rast(nc_path)

first_layer <- nsidc_raster[[1]]

# Create a mask where land (value 2540) is retained, and other values are set to NA
land_mask <- ifel(first_layer == 2540 | first_layer == 2530, 1, NA)  # Set land and coastline to 1, all else to NA
# Plot the land mask to verify
plot(land_mask, col = "brown", main = "Land Mask (Brown)")

# Define the output file path for the land mask raster
# land_mask_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/Land_Mask/Land_Mask.tif"
land_mask_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/Datasets/Land_Mask/land_mask.tif"

# Save the land mask as a raster
writeRaster(land_mask, land_mask_path, overwrite = TRUE)

# Optionally, convert the land mask to polygons (shapefile)
land_mask_shapefile_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/Datasets/Land_Mask/land_mask.shp"

land_polygon <- as.polygons(land_mask, dissolve=TRUE)  # Convert the mask to polygons
writeVector(land_polygon, land_mask_shapefile_path, overwrite = TRUE)

# Load and verify the shapefile if needed
loaded_land_polygon <- vect(land_mask_shapefile_path)
plot(loaded_land_polygon, col = "brown", main = "Land Mask (Shapefile)")
