library(terra)


# Load the raster
raster_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/rast/IBCSO_v2_bed_WGS84_Clip.tif"
rast <- rast(raster_path)


# Define the reclassification matrix for the continental shelf
# The format is c(from, to, newvalue)
reclass_matrix_shelf <- matrix(c(-500, 0, 1,  # Identify the continental shelf
                                 -Inf, -501, NA,  # Exclude deeper than -500m
                                 1, Inf, NA),    # Exclude above sea level
                               byrow = TRUE, ncol = 3)

# Reclassify the raster using the defined matrix
continental_shelf_rast <- classify(rast, reclass_matrix_shelf)

# Specify the output path for the continental shelf
output_path_shelf <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/rast/IBCSO_v2_bed_WGS84_ContinentalShelf.tif"

# Save the raster
writeRaster(continental_shelf_rast, output_path_shelf, overwrite=TRUE)
