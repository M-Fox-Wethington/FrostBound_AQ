# Overview:
# This script loads a series of raster files representing daily sea ice concentration data,
# extracts the date from each filename, assigns the date as the name for each raster layer,
# combines these layers into a single SpatRaster object, and finally demonstrates how to
# convert a raster layer to a dataframe for further analysis with dplyr.

# Load necessary libraries
library(terra) # For working with spatial data, specifically raster data
library(dplyr) # For data manipulation and analysis

# Function to extract the date from the filename using regular expression
# The first implementation of extractDate is redundant and is replaced by the second one.
# extractDate <- function(filename) {
#   matches <- regexec("_(\\d{8})\\.", filename)
#   match <- regmatches(filename, matches)
#   
#   if (length(match[[1]]) > 1) {
#     return(match[[1]][2])
#   } else {
#     return(NA)
#   }
# }

# Define the directory containing the raster files
raster_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSRE-Unified_Winter_Processed_Subset"
# List all .tif raster files in the directory
raster_files <- list.files(path = raster_dir, pattern = "\\.tif$", full.names = TRUE)

# Updated function to extract the date from filenames
# Assumes the date is in YYYYMMDD format right before the file extension
extractDate <- function(filename) {
  gsub(".*_(\\d{8}).tif$", "\\1", basename(filename))
}

# Load raster files and assign names based on the extracted date
rasters <- lapply(raster_files, function(f) {
  r <- rast(f) # Load the raster file
  names(r) <- extractDate(f) # Set the layer name to the extracted date
  return(r)
})

# Combine the individual rasters into a single SpatRaster object for easier handling
raster_stack <- rast(rasters)

# Print the names of the layers in the raster stack to verify correct naming
print(names(raster_stack))

# Example of converting a raster layer to a dataframe for further analysis
# Here, we take the first layer as an example. This conversion facilitates
# the use of dplyr for data manipulation and analysis.
df <- as.data.frame(raster_stack, xy = TRUE, na.rm = TRUE) %>%
  mutate(Date = names(raster_stack)[1]) # Add a 'Date' column based on the layer names


# Define the path for saving the raster stack
output_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSRE-Unified_Winter_Processed_Subset/stack/AMSR-E-Unified_Winter_SIC_Stack_Subset.tif"

# Use writeRaster to save the raster stack
writeRaster(raster_stack, filename = output_path, overwrite = TRUE)

# Load the raster stack back into R
# loaded_raster_stack <- rast(output_path)

# Print summary to check
print(loaded_raster_stack)

str(df)
head(df)

gc()
