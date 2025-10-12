# Load necessary libraries
library(terra)
library(dplyr)

# Load your Gentoo penguin dataset
gentoo_file_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/GentooCounts_48_1.csv"
gentoo_data <- read.csv(gentoo_file_path)
str(gentoo_data)
# Convert Gentoo penguin data to a spatial object
gentoo_spatial <- vect(gentoo_data, geom = c("longitude_epsg_4326", "latitude_epsg_4326"), crs = "EPSG:4326")

# Load the study area shapefile

study_areas <- vect("D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/shp/subregions/Frostbound_AQ_Subregions_EPSG_3976.shp")

# Transform the Gentoo penguin spatial data to the same CRS as the study areas
gentoo_spatial <- project(gentoo_spatial, crs(study_areas))

# Perform the spatial join to append study area information
gentoo_with_study_area <- terra::intersect(gentoo_spatial, study_areas)
output_file_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/GentooCounts_with_StudyArea.shp"
writeVector(gentoo_with_study_area, output_file_path)


# Convert the result back to a data frame
gentoo_with_study_area_df <- as.data.frame(gentoo_with_study_area)

str(gentoo_with_study_area_df)

# Save the updated dataset
output_file_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/GentooCounts_with_StudyArea.csv"
write.csv(gentoo_with_study_area_df, output_file_path, row.names = FALSE)

# View the structure of the updated dataset
str(gentoo_with_study_area_df)

