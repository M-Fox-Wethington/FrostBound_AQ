# Load necessary libraries
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(patchwork)  # For combining plots

# Define file paths for the raster datasets and the coastline mask
amsr_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/amsr_common_dates.tif"
nsidc_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_common_dates.tif"
nsidc_harmonized_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_harmonized_common_dates.tif"
coastline_mask_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/NSIDC_Sea-Ice-Index_Coastline_Mask/Sea-Ice-Index_Study_Area_Coastline-Mask_12km.tif"

# Load the raster datasets and the coastline mask
amsr_stack <- rast(amsr_path)
nsidc_stack <- rast(nsidc_path)
nsidc_harmonized_stack <- rast(nsidc_harmonized_path)
coastline_mask <- rast(coastline_mask_path)

# Extract the time information from the rasters
print("Extracting time information...")
amsr_dates <- time(amsr_stack)
nsidc_dates <- time(nsidc_stack)

# Ensure dates are in Date format and find common dates
print("Finding common dates...")
amsr_dates <- as.Date(amsr_dates)
nsidc_dates <- as.Date(nsidc_dates)
common_dates <- as.Date(intersect(amsr_dates, nsidc_dates))

# Find the indices of common dates in the AMSR and NSIDC datasets
print("Finding indices of common dates...")
amsr_indices <- which(amsr_dates %in% common_dates)
nsidc_indices <- which(nsidc_dates %in% common_dates)

# Subset the layers for the overlapping period
print("Subsetting layers for the overlapping period...")
amsr_overlap <- subset(amsr_stack, amsr_indices)
nsidc_overlap <- subset(nsidc_stack, nsidc_indices)

# Set the layer names to the corresponding common dates
print("Setting layer names for AMSR and NSIDC overlap layers...")
names(amsr_overlap) <- common_dates
names(nsidc_overlap) <- common_dates

# 1. Check alignment: Compare extents, resolution, and CRS, and resample if needed
if (!compareGeom(amsr_overlap, nsidc_overlap, stopOnError = FALSE)) {
  nsidc_overlap <- resample(nsidc_overlap, amsr_overlap)
}
if (!compareGeom(amsr_overlap, coastline_mask, stopOnError = FALSE)) {
  coastline_mask <- resample(coastline_mask, amsr_overlap)
}

# Apply the coastline mask using the mask() function from terra
amsr_masked <- mask(amsr_overlap, coastline_mask, maskvalue = 1)
nsidc_masked <- mask(nsidc_overlap, coastline_mask, maskvalue = 1)

# Continue with the rest of the analysis using the masked datasets
# Calculate the area of one cell in square kilometers
cell_area_sq_km <- prod(res(amsr_masked)) / 1e6  # Adjust this for different rasters if resolutions differ

# Function to calculate sea ice extent statistics
calculate_layer_stats <- function(layer, date, cell_area_sq_km) {
  valid_ice_cells <- terra::global(layer, fun = function(x) sum(x >= 0.15, na.rm = TRUE))  # Explicit call to terra::global
  total_ice_area_sq_km <- (valid_ice_cells * cell_area_sq_km)  # Total ice area in square kilometers
  return(data.frame(ice_extent_km = total_ice_area_sq_km, date = date))
}

# Function to calculate the average pixel-level difference
calculate_pixel_level_diff <- function(layer1, layer2, cell_area_sq_km) {
  abs_diff <- abs(layer1 - layer2)
  pixel_diff_km <- terra::global(abs_diff, fun = function(x) mean(x, na.rm = TRUE)) * cell_area_sq_km  # Adjust for km²
  return(pixel_diff_km)
}

# Function to process the raster data and calculate statistics
process_data <- function(raster, dates, cell_area_sq_km) {
  sea_ice_stats <- lapply(1:nlyr(raster), function(i) {
    calculate_layer_stats(raster[[i]], dates[i], cell_area_sq_km)
  })
  sea_ice_stats_df <- do.call(rbind, sea_ice_stats) %>%
    mutate(Date = as.Date(date), Year = year(Date), Month = month(Date, label = TRUE, abbr = TRUE))
  return(sea_ice_stats_df)
}

# Process the AMSR and NSIDC data (ensure file paths for saving results)
amsr_stats_save <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/amsr_stats_overlap.rds"
nsidc_stats_save <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/nsidc_stats_overlap.rds"

if (file.exists(amsr_stats_save)) {
  amsr_stats <- readRDS(amsr_stats_save)
} else {
  amsr_stats <- process_data(amsr_masked, common_dates, cell_area_sq_km) %>%
    mutate(Source = "AMSR")
  saveRDS(amsr_stats, amsr_stats_save)
}

if (file.exists(nsidc_stats_save)) {
  nsidc_stats <- readRDS(nsidc_stats_save)
} else {
  nsidc_stats <- process_data(nsidc_masked, common_dates, cell_area_sq_km) %>%
    mutate(Source = "NSIDC")
  saveRDS(nsidc_stats, nsidc_stats_save)
}

# Rename 'global' to 'ice_extent_km' in the stats
amsr_stats <- amsr_stats %>% rename(ice_extent_km = global)
nsidc_stats <- nsidc_stats %>% rename(ice_extent_km = global)

# Calculate the pixel-level average differences
pixel_diff <- lapply(1:nlyr(amsr_masked), function(i) {
  calculate_pixel_level_diff(amsr_masked[[i]], nsidc_masked[[i]], cell_area_sq_km)
})

# Unlist and bind the pixel differences to the statistics
pixel_diff <- unlist(pixel_diff)

# Save or load the differences
diff_save <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/differences_overlap.rds"
if (file.exists(diff_save)) {
  differences <- readRDS(diff_save)
} else {
  differences <- bind_rows(amsr_stats, nsidc_stats) %>%
    pivot_wider(names_from = Source, values_from = ice_extent_km) %>%
    mutate(Diff_AMSR_NS = AMSR - NSIDC,
           Avg_Pixel_Diff = pixel_diff)
  saveRDS(differences, diff_save)
}

# Clean up any NAs, NaNs, or Infinities in the differences
differences_filtered <- differences %>%
  filter(!is.na(Diff_AMSR_NS), !is.nan(Diff_AMSR_NS), !is.infinite(Diff_AMSR_NS))

# Summarize the differences
summary_diff <- differences_filtered %>%
  summarise(
    Mean_Diff_AMSR_NS = mean(Diff_AMSR_NS, na.rm = TRUE),
    SD_Diff_AMSR_NS = sd(Diff_AMSR_NS, na.rm = TRUE),
    Mean_Pixel_Diff = mean(Avg_Pixel_Diff, na.rm = TRUE)
  )

# Print summary
print(summary_diff)

# Plot the differences
p_diff <- ggplot(differences_filtered, aes(x = Year, y = Diff_AMSR_NS)) +
  geom_line(color = "#377eb8") +
  labs(title = "Difference in Sea Ice Extent: AMSR vs NSIDC", y = "Difference (sq km)", x = "Year")

# Print the plot
print(p_diff)
