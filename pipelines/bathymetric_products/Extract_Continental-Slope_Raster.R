library(terra)

# Load the raster
raster_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/rast/IBCSO_v2_bed_WGS84_Clip.tif"
rast <- rast(raster_path)

# Define the reclassification matrix
# The format is c(from, to, newvalue)
reclass_matrix <- matrix(c(-Inf, -3001, NA,
                           -3000, -500, 1,
                           -499, Inf, NA), byrow = TRUE, ncol = 3)

# Reclassify the raster
reclassified_rast <- classify(rast, reclass_matrix)

# Specify the output path
output_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/rast/IBCSO_v2_bed_WGS84_Continental-Slope.tif"

# Save the raster
writeRaster(reclassified_rast, output_path, overwrite=TRUE)
