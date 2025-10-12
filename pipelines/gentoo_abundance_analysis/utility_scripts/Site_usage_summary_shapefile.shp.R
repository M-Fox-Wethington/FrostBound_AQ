# Load necessary libraries
library(sf)
library(dplyr)

# Assuming penguin_abundance_filtered is your dataframe and has columns site_id, latitude, longitude, and consecutive_used
# Summarize the number of years each site was used in the model
site_usage_summary <- penguin_abundance_filtered %>%
  filter(consecutive_used == "Used") %>%
  group_by(site_id) %>%
  summarize(
    num_years_used = n(),
    latitude = first(latitude),
    longitude = first(longitude)
  )

# Create an sf object from the summarized data
site_usage_sf <- st_as_sf(
  site_usage_summary,
  coords = c("longitude", "latitude"),
  crs = 4326
)

# Define the output path for the shapefile
output_shapefile_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/gentoo-sites-used/site_usage_summary.shp"

# Write the sf object to a shapefile
st_write(site_usage_sf, output_shapefile_path)

# Print a message indicating the shapefile has been created
cat("Shapefile created at:", output_shapefile_path, "\n")
