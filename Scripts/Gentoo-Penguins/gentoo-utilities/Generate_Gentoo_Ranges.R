# Load the necessary libraries
library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# Load Gentoo penguin data
penguin_data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/modeled_gentoo_parameters.csv")

# Filter to unique site_ids
unique_sites <- penguin_data %>%
  distinct(site_id, .keep_all = TRUE)

# Load coastline shapefile
coastline_path <- "D:/gis_data_antarctica/coastlines/polygon/add_coastline_medium_res_polygon_v7_5.shp"
coastline <- st_read(coastline_path)

# Load study area shapefile
study_area_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/shp/subregions/Frostbound_AQ_Subregions_EPSG_3976.shp"
study_area <- st_read(study_area_path)

# Convert Gentoo penguin data to sf object
gentoo_sf <- st_as_sf(unique_sites, coords = c("longitude", "latitude"), crs = 4326)
ant_proj <- st_crs(study_area)$wkt
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
land_water_rast <- mask(land_water_rast, vect(study_area))

# List to store buffers
buffer_list <- list()

# Loop through unique colonies to calculate buffers
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
  
  # Add site_id to each buffer polygon
  buffer_sf$site_id <- gentoo_sf$site_id[i]
  
  
  buffer_list[[i]] <- buffer_sf
}

# Combine the separate buffers into a single sf object for plotting and saving
combined_buffers_sf <- do.call(rbind, buffer_list)
# 
# # Plot the combined buffers over the study area with enhanced visibility
# ggplot() +
#   geom_sf(data = study_area, fill = NA, color = "black") +
#   geom_sf(data = coastline, fill = "grey") +
#   geom_sf(data = combined_buffers_sf, aes(fill = site_id), alpha = 0.5) +
#   geom_sf(data = gentoo_sf, color = "red", size = 3) +
#   theme_minimal() +
#   labs(title = "Potential Range for All Gentoo Penguin Colonies (300 km)",
#        fill = "Colony") +
#   theme(plot.title = element_text(hjust = 0.5))

# Save the combined buffers to a shapefile
st_write(combined_buffers_sf, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/gentoo_home-ranges/GEPG_homerange_300km.shp", append=FALSE)
