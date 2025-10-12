# Load necessary libraries
library(terra)

# Load the NetCDF file or raster file
nc_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/25km_Sea-Ice-Index/processed_clipped/S_20220307_concentration_v3.0.tif"
nsidc_raster <- rast(nc_path)

first_layer <- nsidc_raster[[1]]  # Select the first layer

# Create a mask where all values are set to 1 (turn the entire raster into 1s)
black_mask <- first_layer * 0 + 1  # Set all cells to 1

# Plot the black mask to verify
plot(black_mask, col = "black", main = "Black Mask (All values set to 1)")

# Define the output file path for the black mask raster
black_mask_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets//Land_Mask/Black_Mask.tif"

# Save the black mask as a raster
writeRaster(black_mask, black_mask_path, overwrite = TRUE)
