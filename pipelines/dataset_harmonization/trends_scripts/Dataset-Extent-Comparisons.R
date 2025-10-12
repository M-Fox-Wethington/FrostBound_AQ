
# Load necessary libraries
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(patchwork)  # For combining plots

# Define file paths for the raster datasets
amsr_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/amsr_common_dates.tif"
nsidc_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_common_dates.tif"
nsidc_harmonized_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_harmonized_common_dates.tif"

# Define file paths for intermediate saves
amsr_stats_save <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/amsr_stats.rds"
nsidc_stats_save <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/nsidc_stats.rds"
nsidc_harmonized_stats_save <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/nsidc_harmonized_stats.rds"
diff_save <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/differences.rds"

# Load the raster datasets
amsr_stack <- rast(amsr_path)
nsidc_stack <- rast(nsidc_path)
nsidc_harmonized_stack <- rast(nsidc_harmonized_path)

# Calculate the area of one cell in square kilometers
cell_area_sq_km <- prod(res(amsr_stack)) / 1e6  # Adjust this for different rasters if resolutions differ

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

# Process AMSR, NSIDC, and Harmonized data (ensure file paths for saving results)
if (file.exists(amsr_stats_save)) {
  amsr_stats <- readRDS(amsr_stats_save)
} else {
  amsr_stats <- process_data(amsr_stack, amsr_dates, cell_area_sq_km) %>%
    mutate(Source = "AMSR")
  saveRDS(amsr_stats, amsr_stats_save)
}

if (file.exists(nsidc_stats_save)) {
  nsidc_stats <- readRDS(nsidc_stats_save)
} else {
  nsidc_stats <- process_data(nsidc_stack, nsidc_dates, cell_area_sq_km) %>%
    mutate(Source = "NSIDC")
  saveRDS(nsidc_stats, nsidc_stats_save)
}

if (file.exists(nsidc_harmonized_stats_save)) {
  nsidc_harmonized_stats <- readRDS(nsidc_harmonized_stats_save)
} else {
  nsidc_harmonized_stats <- process_data(nsidc_harmonized_stack, nsidc_harmonized_dates, cell_area_sq_km) %>%
    mutate(Source = "NSIDC_Harmonized")
  saveRDS(nsidc_harmonized_stats, nsidc_harmonized_stats_save)
}

# Rename 'global' to 'ice_extent_km' in the stats
amsr_stats <- amsr_stats %>% rename(ice_extent_km = global)
nsidc_stats <- nsidc_stats %>% rename(ice_extent_km = global)
nsidc_harmonized_stats <- nsidc_harmonized_stats %>% rename(ice_extent_km = global)

# Calculate the pixel-level average differences
pixel_diff_amsr <- lapply(1:nlyr(nsidc_harmonized_stack), function(i) {
  calculate_pixel_level_diff(nsidc_harmonized_stack[[i]], amsr_stack[[i]], cell_area_sq_km)
})

pixel_diff_nsidc <- lapply(1:nlyr(nsidc_harmonized_stack), function(i) {
  calculate_pixel_level_diff(nsidc_harmonized_stack[[i]], nsidc_stack[[i]], cell_area_sq_km)
})

# Unlist and bind the pixel differences to the statistics
pixel_diff_amsr <- unlist(pixel_diff_amsr)
pixel_diff_nsidc <- unlist(pixel_diff_nsidc)

# Save or load the differences
if (file.exists(diff_save)) {
  differences <- readRDS(diff_save)
} else {
  differences <- bind_rows(amsr_stats, nsidc_stats, nsidc_harmonized_stats) %>%
    pivot_wider(names_from = Source, values_from = ice_extent_km) %>%
    mutate(
      Diff_AMSR = NSIDC_Harmonized - AMSR,
      Diff_NSIDC = NSIDC_Harmonized - NSIDC,
      Avg_Pixel_Diff_AMSR = pixel_diff_amsr,  
      Avg_Pixel_Diff_NSIDC = pixel_diff_nsidc
    )
  saveRDS(differences, diff_save)
}

# Clean up any NAs, NaNs, or Infinities in the differences
differences_filtered <- differences %>%
  filter(!is.na(Diff_AMSR), !is.nan(Diff_AMSR), !is.infinite(Diff_AMSR),
         !is.na(Diff_NSIDC), !is.nan(Diff_NSIDC), !is.infinite(Diff_NSIDC))

# Summarize the differences
summary_diff <- differences_filtered %>%
  summarise(
    Mean_Diff_AMSR = mean(Diff_AMSR, na.rm = TRUE),
    SD_Diff_AMSR = sd(Diff_AMSR, na.rm = TRUE),
    Mean_Diff_NSIDC = mean(Diff_NSIDC, na.rm = TRUE),
    SD_Diff_NSIDC = sd(Diff_NSIDC, na.rm = TRUE),
    Mean_Pixel_Diff_AMSR = mean(Avg_Pixel_Diff_AMSR, na.rm = TRUE),
    Mean_Pixel_Diff_NSIDC = mean(Avg_Pixel_Diff_NSIDC, na.rm = TRUE)
  )

# Print summary
print(summary_diff)

# Plot the differences
p_diff_amsr <- ggplot(differences_filtered, aes(x = Year, y = Diff_AMSR)) +
  geom_line(color = "#377eb8") +
  labs(title = "Difference in Sea Ice Extent: NSIDC Harmonized vs AMSR", y = "Difference (sq km)", x = "Year")

p_diff_nsidc <- ggplot(differences_filtered, aes(x = Year, y = Diff_NSIDC)) +
  geom_line(color = "#ff7f00") +
  labs(title = "Difference in Sea Ice Extent: NSIDC Harmonized vs NSIDC", y = "Difference (sq km)", x = "Year")

# Combine the plots
combined_diff_plot <- p_diff_amsr / p_diff_nsidc + plot_layout(guides = 'collect') + theme(legend.position = "right")

# Print the combined plot
print(combined_diff_plot)
