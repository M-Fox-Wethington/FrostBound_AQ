# Load necessary libraries
library(terra)
library(sf)
library(ggplot2)

# Load Gentoo penguin dataset
gentoo_file_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/GentooCounts_48_1.csv"
gentoo_data <- read.csv(gentoo_file_path)

# Load Antarctic coastline shapefile
coastline_path <- "D:/gis_data_antarctica/coastlines/polygon/add_coastline_medium_res_polygon_v7_5.shp"
coastline <- st_read(coastline_path)

# Load study area shapefile
study_area_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/shp/subregions/Frostbound_AQ_Subregions_EPSG_3976.shp"
study_area <- st_read(study_area_path)

# Convert Gentoo penguin data to sf object
gentoo_sf <- st_as_sf(gentoo_data, coords = c("longitude_epsg_4326", "latitude_epsg_4326"), crs = 4326)

# Transform coordinates to a suitable Antarctic projection
ant_proj <- st_crs(study_area)$wkt  # Use the same projection as the study area
gentoo_sf <- st_transform(gentoo_sf, ant_proj)
coastline <- st_transform(coastline, ant_proj)
study_area <- st_transform(study_area, ant_proj)

# Crop the coastline to the study area
coastline <- st_intersection(coastline, study_area)

# Create a raster with resolution to match the sea ice data (25 km)
resolution <- 25000  # in meters
r <- rast(ext(st_bbox(study_area)), resolution = resolution)
crs(r) <- ant_proj

# Rasterize the coastline (land = NA, water = 0)
coastline_rast <- rasterize(coastline, r, field = 1, background = NA)

# Create a binary land-water raster (land = NA, water = 0)
land_water_rast <- ifel(is.na(coastline_rast), 0, NA)

# Mask the land-water raster with the study area
land_water_rast <- mask(land_water_rast, vect(study_area))

# List to store buffers
buffer_list <- list()

# Loop through all colonies to calculate buffers
for (i in 1:nrow(gentoo_sf)) {
  colony_sf <- gentoo_sf[i, ]
  colony_vect <- vect(colony_sf)
  colony_rast <- rasterize(colony_vect, r, field = 1, background = NA)
  land_water_rast[colony_rast == 1] <- i + 1  # Assign a unique target value to each colony cell
  
  # Calculate distances using gridDist from the target cells
  distance_rast <- gridDist(land_water_rast, target = i + 1)
  
  # Set buffer distance
  buffer_distance <- 300000  # 300 km in meters
  
  # Create a binary raster where cells within 300 km distance are set to 1
  buffer_rast <- distance_rast <= buffer_distance
  
  # Convert the binary buffer raster back to polygons
  buffer_poly <- as.polygons(buffer_rast, dissolve = TRUE)
  
  # Convert to sf object and filter out background layers (layer = 0)
  buffer_sf <- st_as_sf(buffer_poly)
  buffer_sf <- buffer_sf[buffer_sf$layer != 0, ]
  
  buffer_list[[i]] <- buffer_sf
}

# Combine the separate buffers into a single sf object for plotting
combined_buffers_sf <- do.call(rbind, buffer_list)

# Plot the combined buffers over the study area with enhanced visibility
ggplot() +
  geom_sf(data = study_area, fill = NA, color = "black") +
  geom_sf(data = coastline, fill = "grey") +
  geom_sf(data = combined_buffers_sf, aes(fill = as.factor(layer)), alpha = 0.5) +
  geom_sf(data = gentoo_sf, color = "red", size = 3) +
  theme_minimal() +
  labs(title = "Potential Range for All Gentoo Penguin Colonies (300 km)",
       fill = "Colony") +
  theme(plot.title = element_text(hjust = 0.5))
