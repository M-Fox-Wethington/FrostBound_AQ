# Load necessary libraries
library(terra)
library(dplyr)

# Load the AMSR-Unified 12.5 km dataset
amsr_12km <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-Unified_12km/stack/substack/AMSR_12km_Full_Study_Area.nc")

# Load the NSIDC 25 km dataset
nsidc_25km <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/25km_Sea-Ice-Index/stack/substack/NSIDC_25km_Full_Study_Area.nc")

# Verify extents
ext_amsr <- ext(amsr_12km)
ext_nsidc <- ext(nsidc_25km)

print(ext_amsr)
print(ext_nsidc)

# Check origins
origin_amsr <- origin(amsr_12km)
origin_nsidc <- origin(nsidc_25km)

print(origin_amsr)
print(origin_nsidc)

# Align the origin of the NSIDC raster to match the AMSR raster
nsidc_aligned <- shift(nsidc_25km, dx = origin_amsr[1] - origin_nsidc[1], dy = origin_amsr[2] - origin_nsidc[2])

# Verify origins again
origin_nsidc_aligned <- origin(nsidc_aligned)
print(origin_nsidc_aligned)

# Manually define the intersecting extent with ymin rounded down to 174999
intersection_extent <- ext(
  max(ext_amsr$xmin, ext(nsidc_aligned)$xmin),
  min(ext_amsr$xmax, ext(nsidc_aligned)$xmax),
  174999,
  min(ext_amsr$ymax, ext(nsidc_aligned)$ymax)
)

print(intersection_extent)

# Crop both datasets to the intersecting extent
amsr_cropped <- crop(amsr_12km, intersection_extent)
nsidc_cropped <- crop(nsidc_aligned, intersection_extent)

# Ensure the cropped extents are the same by aligning extents precisely
ext_amsr_cropped <- ext(amsr_cropped)
ext_nsidc_cropped <- ext(nsidc_cropped)

print(ext_amsr_cropped)
print(ext_nsidc_cropped)


# Resample the cropped NSIDC dataset to 12.5 km resolution to match the AMSR dataset using bilinear interpolation
nsidc_resampled_to_12_5km <- resample(nsidc_cropped, amsr_cropped, method = "bilinear")



# Normalize raster values to a 0-1 range
normalize_raster <- function(r) {
  r_min <- global(r, min, na.rm = TRUE)[[1]]
  r_max <- global(r, max, na.rm = TRUE)[[1]]
  (r - r_min) / (r_max - r_min)
}

amsr_normalized <- normalize_raster(amsr_cropped)
nsidc_normalized <- normalize_raster(nsidc_resampled_to_12_5km)

# Extract the time information
amsr_dates <- time(amsr_normalized)
nsidc_dates <- time(nsidc_normalized)

# Ensure dates are in Date format and find common dates
amsr_dates <- as.Date(amsr_dates)
nsidc_dates <- as.Date(nsidc_dates)
common_dates <- intersect(amsr_dates, nsidc_dates)

# Convert the common_dates to Date format
common_dates <- as.Date(common_dates, origin = "1970-01-01")

# Find the indices of common dates in the AMSR and NSIDC datasets
amsr_indices <- which(amsr_dates %in% common_dates)
nsidc_indices <- which(nsidc_dates %in% common_dates)

# Subset the layers for the overlapping period
amsr_overlap <- subset(amsr_normalized, amsr_indices)
nsidc_overlap <- subset(nsidc_normalized, nsidc_indices)

# Harmonize the NSIDC layers to AMSR statistics
harmonize_nsidc_layer <- function(layer_index) {
  nsidc_layer <- nsidc_overlap[[layer_index]]
  amsr_layer <- amsr_overlap[[layer_index]]
  
  # Preserve NA values
  na_mask <- is.na(nsidc_layer)
  
  # Calculate the mean and standard deviation for both layers
  mean_nsidc <- global(nsidc_layer, mean, na.rm = TRUE)[[1]]
  sd_nsidc <- global(nsidc_layer, sd, na.rm = TRUE)[[1]]
  mean_amsr <- global(amsr_layer, mean, na.rm = TRUE)[[1]]
  sd_amsr <- global(amsr_layer, sd, na.rm = TRUE)[[1]]
  
  # Harmonize the NSIDC layer to the AMSR statistics
  harmonized_layer <- (nsidc_layer - mean_nsidc) * (sd_amsr / sd_nsidc) + mean_amsr
  harmonized_layer[na_mask] <- NA
  
  # Preserve NA values
  harmonized_layer[is.na(nsidc_layer)] <- NA
  
  # Apply the zero mask
  zero_mask <- (amsr_layer == 0) & (nsidc_layer == 0)
  harmonized_layer[zero_mask] <- 0
  
  # Clamp the values to the 0-1 range
  clamp(harmonized_layer, lower = 0, upper = 1)
  
  # Clamp the values to the 0-1 range
  clamped_layer <- clamp(harmonized_layer, lower = 0, upper = 1)
  
  clamped_layer
}

# Apply the harmonization to all NSIDC layers
harmonized_nsidc_layers <- lapply(1:nlyr(nsidc_overlap), harmonize_nsidc_layer)

# Create a SpatRaster stack from the harmonized layers
harmonized_nsidc_stack <- rast(harmonized_nsidc_layers)

# Assign common dates as names of the layers in the harmonized NSIDC stack
names(harmonized_nsidc_stack) <- as.character(common_dates)


# # Set the time dimension for the harmonized NSIDC stack
# time(harmonized_nsidc_stack) <- common_dates

# Print the extent and resolution of the harmonized NSIDC stack
print(ext(harmonized_nsidc_stack))
print(res(harmonized_nsidc_stack))

# Output the harmonized raster stack to the specified directory
output_filepath <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/linear-transformation/harmonized_nsidc_stack.tif"
writeRaster(harmonized_nsidc_stack, filename = output_filepath, overwrite = TRUE)

# Print confirmation message
print(paste("Harmonized raster stack saved to:", output_filepath))



# Extract the raster layer
raster_layer <- harmonized_nsidc_stack[[10]]

# Plot the raster layer with a default color palette
plot(raster_layer, col = terrain.colors(100), legend = TRUE)

# Create a mask for NA values
na_mask <- is.na(raster_layer)

# Overlay the NA values in red
plot(na_mask, col = c(NA, "red"), add = TRUE, legend = FALSE)

