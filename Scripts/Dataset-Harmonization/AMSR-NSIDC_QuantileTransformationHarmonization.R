# Load necessary libraries
if (!require(terra)) install.packages('terra')
if (!require(dplyr)) install.packages('dplyr')
if (!require(reshape2)) install.packages('reshape2')
if (!require(tidyr)) install.packages('tidyr')

library(terra)
library(dplyr)
library(reshape2)
library(tidyr)

# Function to normalize raster values to a 0-1 range
normalize_raster <- function(r) {
  r_min <- global(r, min, na.rm = TRUE)[[1]]
  r_max <- global(r, max, na.rm = TRUE)[[1]]
  
  if (is.finite(r_min) && is.finite(r_max)) {
    return((r - r_min) / (r_max - r_min))
  } else {
    warning("Raster contains only NA values or min/max are not finite.")
    return(r)
  }
}

# Load the AMSR-Unified 12.5 km dataset
print("Loading AMSR-Unified 12.5 km dataset...")
amsr_12km <- tryCatch({
  rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-Unified_12km/stack/substack/AMSR_12km_Full_Study_Area.nc")
}, error = function(e) {
  stop("Error loading AMSR-Unified dataset: ", e)
})

# Load the NSIDC 25 km dataset (full timeline from 1979 to 2024)
print("Loading NSIDC 25 km dataset...")
nsidc_25km <- tryCatch({
  rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/25km_Sea-Ice-Index/stack/substack/NSIDC_25km_Full_Study_Area.nc")
}, error = function(e) {
  stop("Error loading NSIDC dataset: ", e)
})

# Verify extents
print("Verifying extents...")
ext_amsr <- ext(amsr_12km)
ext_nsidc <- ext(nsidc_25km)

print(ext_amsr)
print(ext_nsidc)

# Check origins
print("Checking origins...")
origin_amsr <- origin(amsr_12km)
origin_nsidc <- origin(nsidc_25km)

print(origin_amsr)
print(origin_nsidc)

# Align the origin of the NSIDC raster to match the AMSR raster
print("Aligning NSIDC raster to AMSR origin...")
nsidc_aligned <- shift(nsidc_25km, dx = origin_amsr[1] - origin_nsidc[1], dy = origin_amsr[2] - origin_nsidc[2])

# Verify origins again
print("Verifying aligned origins...")
origin_nsidc_aligned <- origin(nsidc_aligned)
print(origin_nsidc_aligned)

# Define the intersecting extent with ymin rounded down to 174999
print("Defining intersecting extent...")
intersection_extent <- ext(
  max(ext_amsr$xmin, ext(nsidc_aligned)$xmin),
  min(ext_amsr$xmax, ext(nsidc_aligned)$xmax),
  174999,
  min(ext_amsr$ymax, ext(nsidc_aligned)$ymax)
)

print(intersection_extent)

# Crop both datasets to the intersecting extent
print("Cropping datasets to intersecting extent...")
amsr_cropped <- crop(amsr_12km, intersection_extent)
nsidc_cropped <- crop(nsidc_aligned, intersection_extent)

# Verify cropped extents
print("Verifying cropped extents...")
ext_amsr_cropped <- ext(amsr_cropped)
ext_nsidc_cropped <- ext(nsidc_cropped)

print(ext_amsr_cropped)
print(ext_nsidc_cropped)

# Resample the cropped NSIDC dataset to 12.5 km resolution using bilinear interpolation
print("Resampling NSIDC dataset to 12.5 km resolution...")
nsidc_resampled_to_12_5km <- resample(nsidc_cropped, amsr_cropped, method = "bilinear")

# Normalize each layer individually and then combine
print("Normalizing AMSR raster values...")
amsr_normalized <- rast(lapply(1:nlyr(amsr_cropped), function(i) {
  print(paste("Normalizing AMSR layer", i))
  normalize_raster(amsr_cropped[[i]])
}))

print("Normalizing NSIDC raster values...")
nsidc_normalized <- rast(lapply(1:nlyr(nsidc_resampled_to_12_5km), function(i) {
  print(paste("Normalizing NSIDC layer", i))
  normalize_raster(nsidc_resampled_to_12_5km[[i]])
}))

# Extract the time information
print("Extracting time information...")
amsr_dates <- time(amsr_normalized)
nsidc_dates <- time(nsidc_normalized)

# Ensure dates are in Date format and find common dates
print("Finding common dates...")
amsr_dates <- as.Date(amsr_dates)
nsidc_dates <- as.Date(nsidc_dates)
common_dates <- intersect(amsr_dates, nsidc_dates)

# Find the indices of common dates in the AMSR and NSIDC datasets
print("Finding indices of common dates...")
amsr_indices <- which(amsr_dates %in% common_dates)
nsidc_indices <- which(nsidc_dates %in% common_dates)

# Subset the layers for the overlapping period
print("Subsetting layers for the overlapping period...")
amsr_overlap <- subset(amsr_normalized, amsr_indices)
nsidc_overlap <- subset(nsidc_normalized, nsidc_indices)

# Set the layer names to the corresponding dates
print("Setting layer names for AMSR and NSIDC overlap layers...")
names(amsr_overlap) <- common_dates
names(nsidc_overlap) <- common_dates

# Flatten the rasters to data frames
amsr_df <- as.data.frame(amsr_overlap, xy = TRUE)
nsidc_df <- as.data.frame(nsidc_overlap, xy = TRUE)

# Merge the data frames on the coordinates
merged_df <- merge(amsr_df, nsidc_df, by = c("x", "y"))

# Melt the dataframe to long format
merged_df_long <- melt(merged_df, id.vars = c("x", "y"))

# Separate the date and dataset type from the variable column
merged_df_long <- merged_df_long %>%
  mutate(date = as.Date(gsub("\\..*", "", variable), format = "%Y-%m-%d"),
         dataset = ifelse(grepl("\\.x$", variable), "amsr", "nsidc")) %>%
  select(-variable)

# Spread the dataframe to have separate columns for NSIDC and AMSR
merged_df_wide <- merged_df_long %>%
  pivot_wider(names_from = dataset, values_from = value)

# View the structure of the melted and spread dataframe
str(merged_df_wide)

# Ensure there are no missing values for the regression
merged_df_wide_clean <- merged_df_wide %>%
  filter(!is.na(nsidc) & !is.na(amsr))

# Build the linear regression model using the cleaned data
lm_model <- lm(amsr ~ nsidc, data = merged_df_wide_clean)

# Summary of the linear model
summary(lm_model)

# Function to apply linear model to a raster stack using the app function from terra
apply_lm_to_raster <- function(raster_stack, lm_model) {
  # Preserve the original time attribute
  original_time <- time(raster_stack)
  
  processed_stack <- app(raster_stack, fun = function(x) {
    # Apply the linear model
    predicted_values <- lm_model$coefficients[1] + lm_model$coefficients[2] * x
    # Clamp predicted values to a 0-100 range
    predicted_values <- pmin(pmax(predicted_values, 0), 100)
    # Set all values below 15 to 0
    predicted_values[predicted_values < 15] <- 0
    return(predicted_values)
  })
  
  # Restore the original time attribute
  time(processed_stack) <- original_time
  names(processed_stack) <- names(raster_stack)
  return(processed_stack)
}

# Apply the linear model to the entire NSIDC dataset
processed_stack <- apply_lm_to_raster(nsidc_resampled_to_12_5km, lm_model)

# Extract the first layer for plotting
nsidc_first_layer <- nsidc_resampled_to_12_5km[[4962]]
processed_first_layer <- processed_stack[[4962]]
amsr <- amsr_overlap[[1432]]



# Plot the processed layer to verify
par(mfrow = c(1, 3))  # Set up the plotting area for side-by-side plots

plot(nsidc_first_layer, main = "NSIDC Aligned - First Layer")
plot(processed_first_layer, main = "Processed Layer - First Layer")
plot(amsr, main = "AMSR Layer - First Layer")
# Reset the plotting area
par(mfrow = c(1, 1))

