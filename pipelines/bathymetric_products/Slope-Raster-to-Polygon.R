# Path to the raster
raster_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/bathymetry/continental_slope/rast/IBSCO_V2_Continental-Slope.tif"

# Load the raster
slope_rast <- rast(raster_path)

# Mask out all other values except 1
slope_masked <- slope_rast == 1

# Convert raster to polygons
slope_polygons <- as.polygons(slope_masked, dissolve = TRUE)


plot(slope_polygons)


# Specify the output path for the polygon shapefile
output_poly_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/bathymetry/continental_slope/shp/Continental_Slope_Polygon.shp"

# Save the polygon shapefile
writeVector(slope_polygons, output_poly_path, overwrite=TRUE)
