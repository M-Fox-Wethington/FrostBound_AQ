# Load necessary libraries
library(terra)

# Define file paths
amsr_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/normalized_rasters_2012-2023/amsr_normalized_2012-2023.tif"
nsidc_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/normalized_rasters_2012-2023/nsidc_normalized_2012-2023.tif"
nsidc_harmonized_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/complete-harmonized-dataset/nsidc_12_5km_harmonized_1979-2023.tif"

# Load datasets
amsr_stack <- rast(amsr_path)
nsidc_stack <- rast(nsidc_path)
nsidc_harmonized_stack <- rast(nsidc_harmonized_path)

# Extract the date information from the raster metadata
amsr_dates <- as.Date(terra::time(amsr_stack), origin = "1970-01-01")
nsidc_dates <- as.Date(terra::time(nsidc_stack), origin = "1970-01-01")
nsidc_harmonized_dates <- as.Date(terra::time(nsidc_harmonized_stack), origin = "1970-01-01")

# Find common dates between AMSR and NSIDC Harmonized datasets
common_dates <- as.Date(intersect(amsr_dates, nsidc_harmonized_dates))

# Find indices for common dates in the NSIDC Harmonized stack
nsidc_harmonized_indices <- which(nsidc_harmonized_dates %in% common_dates)

# Subset the NSIDC Harmonized stack for the overlapping period
cat("Subsetting NSIDC Harmonized stack for the overlapping period...\n")
nsidc_harmonized_overlap <- subset(nsidc_harmonized_stack, nsidc_harmonized_indices)

# Find indices for common dates in the NSIDC stack
nsidc_indices <- which(nsidc_dates %in% common_dates)

# Subset the NSIDC stack for the overlapping period
cat("Subsetting NSIDC stack for the overlapping period...\n")
nsidc_overlap <- subset(nsidc_stack, nsidc_indices)

# Define the target directory for exporting
output_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/harmonized-dataset-comparison"

# Export the subsetted rasters
cat("Exporting subsetted rasters...\n")
writeRaster(nsidc_harmonized_overlap, filename = file.path(output_dir, "nsidc_harmonized_overlap.tif"), overwrite = TRUE)
writeRaster(nsidc_overlap, filename = file.path(output_dir, "nsidc_overlap.tif"), overwrite = TRUE)

# Print the structure of each subsetted dataset
print("AMSR 12.5 km Dataset Structure:")
print(list(dimensions = dim(amsr_stack), extent = ext(amsr_stack), crs = crs(amsr_stack), nlyr = nlyr(amsr_stack), resolution = res(amsr_stack)))

print("NSIDC 25 km Overlap Dataset Structure:")
print(list(dimensions = dim(nsidc_overlap), extent = ext(nsidc_overlap), crs = crs(nsidc_overlap), nlyr = nlyr(nsidc_overlap), resolution = res(nsidc_overlap)))

print("NSIDC Harmonized 12.5 km Overlap Dataset Structure:")
print(list(dimensions = dim(nsidc_harmonized_overlap), extent = ext(nsidc_harmonized_overlap), crs = crs(nsidc_harmonized_overlap), nlyr = nlyr(nsidc_harmonized_overlap), resolution = res(nsidc_harmonized_overlap)))
