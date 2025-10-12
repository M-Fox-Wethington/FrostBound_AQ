# Load necessary libraries
library(terra)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)  # For combining plots

# Define file paths for the raster datasets
amsr_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/amsr_common_dates.tif"
nsidc_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_common_dates.tif"
nsidc_harmonized_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_harmonized_common_dates.tif"

# Load the raster datasets
amsr_stack <- rast(amsr_path)
nsidc_stack <- rast(nsidc_path)
nsidc_harmonized_stack <- rast(nsidc_harmonized_path)

# Extract the date information from the raster metadata
amsr_dates <- as.Date(terra::time(amsr_stack), origin = "1970-01-01")
nsidc_dates <- as.Date(terra::time(nsidc_stack), origin = "1970-01-01")
nsidc_harmonized_dates <- as.Date(terra::time(nsidc_harmonized_stack), origin = "1970-01-01")

# Calculate the area of one cell in square kilometers
cell_area_sq_km <- prod(res(amsr_stack)) / 1e6

# Function to calculate sea ice extent statistics
calculate_layer_stats <- function(layer, date, cell_area_sq_km) {
  valid_ice_cells <- global(layer, fun = function(x) sum(x >= 0.15, na.rm = TRUE))
  total_ice_area_sq_km <- (valid_ice_cells * cell_area_sq_km)  # Total ice area in square kilometers
  return(data.frame(ice_extent_km = total_ice_area_sq_km, date = date))
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

# Process each raster stack
amsr_stats <- process_data(amsr_stack, amsr_dates, cell_area_sq_km)
nsidc_stats <- process_data(nsidc_stack, nsidc_dates, cell_area_sq_km)
nsidc_harmonized_stats <- process_data(nsidc_harmonized_stack, nsidc_harmonized_dates, cell_area_sq_km)

# Ensure consistent column names for all datasets before merging
amsr_stats <- process_data(amsr_stack, amsr_dates, cell_area_sq_km) %>%
  mutate(Source = "AMSR") %>%
  rename(ice_extent_km = global)  # Assuming 'global' is the sea ice extent column

nsidc_stats <- process_data(nsidc_stack, nsidc_dates, cell_area_sq_km) %>%
  mutate(Source = "NSIDC") %>%
  rename(ice_extent_km = global)

nsidc_harmonized_stats <- process_data(nsidc_harmonized_stack, nsidc_harmonized_dates, cell_area_sq_km) %>%
  mutate(Source = "NSIDC_Harmonized") %>%
  rename(ice_extent_km = global)



# Calculate annual statistics using the correct data structure
stats_summary <- bind_rows(amsr_stats, nsidc_stats, nsidc_harmonized_stats) %>%
  group_by(Year, Source) %>%
  summarise(
    mean_extent = mean(ice_extent_km, na.rm = TRUE),
    median_extent = median(ice_extent_km, na.rm = TRUE),
    max_extent = max(ice_extent_km, na.rm = TRUE),
    min_extent = min(ice_extent_km, na.rm = TRUE),
    .groups = "drop"  # Drop grouping for further manipulations
  )

# Define color-blind friendly colors
line_colors <- c("AMSR" = "#377eb8", "NSIDC" = "#ff7f00", "NSIDC Harmonized" = "#e41a1c")

# Determine the y-axis range
y_range <- range(stats_summary$mean_extent, stats_summary$median_extent, 
                 stats_summary$max_extent, stats_summary$min_extent, na.rm = TRUE)

# Plotting each statistic with all sources on each plot
p_mean <- ggplot(stats_summary, aes(x = Year, y = mean_extent, color = Source)) +
  geom_line(size = 1) +
  scale_color_manual(values = line_colors) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Mean Annual Sea Ice Extent", y = "Mean Extent (sq km)")

p_median <- ggplot(stats_summary, aes(x = Year, y = median_extent, color = Source)) +
  geom_line(size = 1) +
  scale_color_manual(values = line_colors) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Median Annual Sea Ice Extent", y = "Median Extent (sq km)")

p_max <- ggplot(stats_summary, aes(x = Year, y = max_extent, color = Source)) +
  geom_line(size = 1) +
  scale_color_manual(values = line_colors) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Maximum Annual Sea Ice Extent", y = "Max Extent (sq km)")

p_min <- ggplot(stats_summary, aes(x = Year, y = min_extent, color = Source)) +
  geom_line(size = 1) +
  scale_color_manual(values = line_colors) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Minimum Annual Sea Ice Extent", y = "Min Extent (sq km)")

# Combine the plots into a grid
combined_plot <- (p_mean | p_median) / (p_max | p_min) + 
  plot_layout(guides = 'collect') +
  theme(legend.position = "right")

# Print the combined plot
print(combined_plot)
