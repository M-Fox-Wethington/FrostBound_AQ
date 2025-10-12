# Load necessary library
library(terra)

# Define file paths
amsr_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/normalized_rasters_2012-2023/amsr_normalized_2012-2023.tif"
nsidc_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/normalized_rasters_2012-2023/nsidc_normalized_2012-2023.tif"
nsidc_harmonized_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_12_5km_harmonized_winter_2012-2023.tif"

# Load datasets
amsr_stack <- rast(amsr_path)
nsidc_stack <- rast(nsidc_path)
nsidc_harmonized_stack <- rast(nsidc_harmonized_path)

# Extract the date information from the raster metadata
amsr_dates <- as.Date(terra::time(amsr_stack), origin = "1970-01-01")
nsidc_dates <- as.Date(terra::time(nsidc_stack), origin = "1970-01-01")
nsidc_harmonized_dates <- as.Date(terra::time(nsidc_harmonized_stack), origin = "1970-01-01")

# Find common dates shared among all three datasets
common_dates <- as.Date(Reduce(intersect, list(amsr_dates, nsidc_dates, nsidc_harmonized_dates)))

# Find indices for common dates in each stack
amsr_indices <- which(amsr_dates %in% common_dates)
nsidc_indices <- which(nsidc_dates %in% common_dates)
nsidc_harmonized_indices <- which(nsidc_harmonized_dates %in% common_dates)

# Subset each stack for the overlapping period
amsr_overlap <- subset(amsr_stack, amsr_indices)
nsidc_overlap <- subset(nsidc_stack, nsidc_indices)
nsidc_harmonized_overlap <- subset(nsidc_harmonized_stack, nsidc_harmonized_indices)


# Define output paths
output_path_amsr <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/amsr_common_dates.tif"
output_path_nsidc <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_common_dates.tif"
output_path_nsidc_harmonized <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison/nsidc_harmonized_common_dates.tif"

# Export the updated raster stacks to the specified directories
writeRaster(amsr_overlap, output_path_amsr, overwrite = TRUE)
writeRaster(nsidc_overlap, output_path_nsidc, overwrite = TRUE)
writeRaster(nsidc_harmonized_overlap, output_path_nsidc_harmonized, overwrite = TRUE)

cat("All common date subsets created and saved successfully.\n")