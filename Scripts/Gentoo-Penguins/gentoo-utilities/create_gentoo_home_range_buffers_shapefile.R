# Load the necessary libraries
library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# Define the function
create_gentoo_buffers <- function(penguin_data_path, sites_used_path, coastline_path, study_area_path, output_shapefile_path, cell_resolution = 25000, buffer_distances = c(50, 100, 150, 250, 500), dissolve_buffers = TRUE) {
  
  # Load Gentoo penguin data
  penguin_data <- read.csv(penguin_data_path)
  
  # Load sites used in the study
  sites_used <- st_read(sites_used_path)
  
  # Filter to unique site_ids that are included in the study
  unique_sites <- penguin_data %>%
    filter(site_id %in% sites_used$site_id) %>%
    distinct(site_id, .keep_all = TRUE)
  
  # Load coastline shapefile
  coastline <- st_read(coastline_path)
  
  # Load study area shapefile
  study_area <- st_read(study_area_path)
  
  # Convert Gentoo penguin data to sf object
  gentoo_sf <- st_as_sf(unique_sites, coords = c("longitude", "latitude"), crs = 4326)
  ant_proj <- st_crs(study_area)$wkt
  gentoo_sf <- st_transform(gentoo_sf, ant_proj)
  coastline <- st_transform(coastline, ant_proj)
  study_area <- st_transform(study_area, ant_proj)
  
  # Extend the bounding box slightly beyond the extent of the coastline
  bbox <- st_bbox(coastline)
  expand_factor <- 0.1
  bbox_expanded <- bbox + c(-1, -1, 1, 1) * expand_factor * (bbox[3:4] - bbox[1:2])
  
  # Create a raster with resolution to match the sea ice data (default 25 km)
  r <- rast(ext(bbox_expanded), resolution = cell_resolution)
  crs(r) <- ant_proj
  
  # Rasterize the coastline (land = NA, water = 0)
  coastline_rast <- rasterize(coastline, r, field = 1, background = NA)
  
  # Create a binary land-water raster (land = NA, water = 0)
  land_water_rast <- ifel(is.na(coastline_rast), 0, NA)
  
  for (buffer_distance in buffer_distances) {
    # List to store buffers
    buffer_list <- list()
    
    # Loop through unique colonies to calculate buffers
    for (i in 1:nrow(gentoo_sf)) {
      print(paste("Processing colony:", gentoo_sf$site_id[i]))  # Debug print
      colony_sf <- gentoo_sf[i, ]
      colony_vect <- vect(colony_sf)
      colony_rast <- rasterize(colony_vect, r, field = 1, background = NA)
      land_water_rast[colony_rast == 1] <- i + 1  # Assign a unique target value to each colony cell
      
      # Calculate distances using gridDist from the target cells
      distance_rast <- gridDist(land_water_rast, target = i + 1)
      
      # Create a binary raster where cells within buffer_distance km are set to 1
      buffer_rast <- distance_rast <= (buffer_distance * 1000)
      
      # Convert the binary buffer raster back to polygons
      buffer_poly <- as.polygons(buffer_rast, dissolve = TRUE)
      
      # Convert to sf object and filter out background layers (layer = 0)
      buffer_sf <- st_as_sf(buffer_poly)
      buffer_sf <- buffer_sf[buffer_sf$layer != 0, ]
      
      # Check if buffer_sf has valid geometries
      if (nrow(buffer_sf) > 0) {
        # Add site_id to each buffer polygon
        buffer_sf$site_id <- gentoo_sf$site_id[i]
        buffer_list[[i]] <- buffer_sf
      } else {
        print(paste("No valid geometry for colony:", gentoo_sf$site_id[i]))  # Debug print
      }
    }
    
    # Combine the separate buffers into a single sf object
    combined_buffers_sf <- do.call(rbind, buffer_list)
    
    if (dissolve_buffers) {
      combined_buffers_sf <- st_union(combined_buffers_sf)
    }
    
    # Save the combined buffers to a shapefile
    output_path <- paste0(output_shapefile_path, "_", buffer_distance, "km.shp")
    st_write(combined_buffers_sf, output_path, append = FALSE)
  }
  
  print("Buffers created and saved.")
}

# Example call to the function
create_gentoo_buffers(
  penguin_data_path = "D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/modeled_gentoo_parameters.csv",
  sites_used_path = "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/gentoo-sites-used/site_usage_summary.shp",
  coastline_path = "D:/gis_data_antarctica/coastlines/polygon/add_coastline_medium_res_polygon_v7_5.shp",
  study_area_path = "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/shp/subregions/Frostbound_AQ_Subregions_EPSG_3976.shp",
  output_shapefile_path = "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/gentoo-home-ranges-updated/gepe_home_ranges",
  cell_resolution = 1000,
  buffer_distances = c(25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500),
  dissolve_buffers = FALSE
)
