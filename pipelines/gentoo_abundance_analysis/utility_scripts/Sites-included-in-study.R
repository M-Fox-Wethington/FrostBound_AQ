# Load the necessary libraries
library(terra)
library(sf)
library(dplyr)

# Load Gentoo penguin data
penguin_data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/gentoo_presence_absence_assumptions.csv")
penguin_abundance_data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/modeled_gentoo_parameters.csv")

# Adjust the season column in the penguin abundance data
penguin_abundance_data <- penguin_abundance_data %>%
  mutate(year = 1970 + season - 1)

# Filter penguin data to keep only consecutively censused sites and clean duplicates
consecutively_censused <- penguin_data %>%
  arrange(site_id, season) %>%
  group_by(site_id) %>%
  mutate(next_season = lead(season),
         next_presence = lead(presence)) %>%
  filter(presence == 1 & next_presence == 1 & season != next_season) %>%
  ungroup() %>%
  select(site_id, season, presence, next_season)

# Merge with penguin abundance data
penguin_abundance_filtered <- penguin_abundance_data %>%
  left_join(consecutively_censused, by = c("site_id", "year" = "season")) %>%
  filter(!is.na(next_season))

# Load study area shapefile
study_area_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/shp/subregions/Frostbound_AQ_Subregions_EPSG_3976.shp"
study_area <- st_read(study_area_path)

# Convert the filtered penguin data to a spatial dataframe
penguin_sites <- penguin_abundance_filtered %>%
  select(site_id, site_name, latitude, longitude) %>%
  distinct() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Reproject the spatial dataframe to match the study area CRS
penguin_sites <- st_transform(penguin_sites, crs = st_crs(study_area))

# Save the spatial dataframe as a shapefile
output_path <- "D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/used_sites.shp"
st_write(penguin_sites, output_path)

# Print message
print(paste("Shapefile created at:", output_path))
