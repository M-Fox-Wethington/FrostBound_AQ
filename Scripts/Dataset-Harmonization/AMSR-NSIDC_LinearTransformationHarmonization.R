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

# Define the manual overlapping extent
manual_overlap_extent <- ext(
  max(ext_amsr$xmin, ext_nsidc$xmin),
  min(ext_amsr$xmax, ext_nsidc$xmax),
  max(ext_amsr$ymin, ext_nsidc$ymin),
  min(ext_amsr$ymax, ext_nsidc$ymax)
)

# Crop both datasets to the manual overlapping extent
amsr_cropped <- crop(amsr_12km, manual_overlap_extent)
nsidc_cropped <- crop(nsidc_25km, manual_overlap_extent)

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
  
  # Clamp the values to the 0-1 range
  clamped_layer <- clamp(harmonized_layer, lower = 0, upper = 1)
  
  clamped_layer
}

# Apply the harmonization to all NSIDC layers
harmonized_nsidc_layers <- lapply(1:nlyr(nsidc_overlap), harmonize_nsidc_layer)

# Create a SpatRaster stack from the harmonized layers
harmonized_nsidc_stack <- rast(harmonized_nsidc_layers)

# Set the time dimension for the harmonized NSIDC stack
time(harmonized_nsidc_stack) <- common_dates

# Print the extent and resolution of the harmonized NSIDC stack
print(ext(harmonized_nsidc_stack))
print(res(harmonized_nsidc_stack))

# Plot a random sample of layers from the original and harmonized datasets
set.seed(123) # For reproducibility
sample_indices <- sample(1:length(common_dates), 10)

plot_layers_side_by_side <- function(index) {
  par(mfrow = c(3, 1))
  
  # Plot AMSR layer and its summary
  hist(values(amsr_overlap[[index]]), main = paste("Histogram of AMSR Layer for", common_dates[index]), xlab = "Value", breaks = 50, col = "blue")
  summary_amsr <- summary(amsr_overlap[[index]])
  print(paste("Summary of AMSR Layer for", common_dates[index]))
  print(summary_amsr)
  
  # Plot NSIDC layer and its summary
  hist(values(nsidc_overlap[[index]]), main = paste("Histogram of NSIDC Layer for", common_dates[index]), xlab = "Value", breaks = 50, col = "green")
  summary_nsidc <- summary(nsidc_overlap[[index]])
  print(paste("Summary of NSIDC Layer for", common_dates[index]))
  print(summary_nsidc)
  
  # Plot Harmonized NSIDC layer and its summary
  hist(values(harmonized_nsidc_stack[[index]]), main = paste("Histogram of Harmonized NSIDC Layer for", common_dates[index]), xlab = "Value", breaks = 50, col = "red")
  summary_harmonized <- summary(harmonized_nsidc_stack[[index]])
  print(paste("Summary of Harmonized NSIDC Layer for", common_dates[index]))
  print(summary_harmonized)
  
  # Plot the layers side by side for comparison
  par(mfrow = c(1, 3))
  plot(amsr_overlap[[index]], main = paste("AMSR Layer for", common_dates[index]))
  plot(nsidc_overlap[[index]], main = paste("NSIDC Layer for", common_dates[index]))
  plot(harmonized_nsidc_stack[[index]], main = paste("Harmonized NSIDC Layer for", common_dates[index]))
  plot(is.na(harmonized_nsidc_stack[[index]]), col = c(NA, "red"), add = TRUE, legend = FALSE)
}

for (i in sample_indices) {
  plot_layers_side_by_side(i)
}
