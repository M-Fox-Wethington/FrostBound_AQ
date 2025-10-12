#THIS IS PERSISTENCE




# Load necessary libraries
library(terra)
library(dplyr)
library(lubridate)

# Define the directory containing the subregion files
chunk_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/25km_Sea-Ice-Index/stack/substack"
chunk_files <- list.files(chunk_dir, pattern = "\\.nc$", full.names = TRUE)

# Define the threshold for sea ice concentration
ice_threshold <- 15

# Define the date range for analysis
start_date <- as.Date("2021-01-01")
end_date <- as.Date("2022-12-31")

# Simplified function to process each chunk file using app
process_chunk <- function(file_path) {
  region_name <- gsub(".*/|\\.nc$", "", file_path)  # Extract the region name from the file path
  message(paste("Processing file:", file_path))
  
  # Load the NSIDC sea ice concentration data for the subregion
  nsidc <- rast(file_path)
  print("Loaded NSIDC data:")
  print(nsidc)
  
  # Filter the sea ice data from 2021 to 2022
  nsidc <- subset(nsidc, which(time(nsidc) >= start_date & time(nsidc) <= end_date))
  print("Filtered NSIDC data:")
  print(nsidc)
  
  # Count the number of days with ice concentration above the threshold for each pixel
  valid_obs <- app(nsidc, function(x) sum(x > ice_threshold, na.rm = TRUE))
  names(valid_obs) <- "Count_Above_Threshold"
  
  # Calculate the statistics for each annual raster layer
  stats <- data.frame(
    Region = region_name,
    Min_Days = min(valid_obs[], na.rm = TRUE),
    Max_Days = max(valid_obs[], na.rm = TRUE),
    Mean_Days = mean(valid_obs[], na.rm = TRUE),
    Median_Days = median(valid_obs[], na.rm = TRUE),
    Start_Date = start_date,
    End_Date = end_date
  )
  
  # Plot the raster
  plot(valid_obs, main = paste("Sea Ice Concentration Above Threshold for", region_name), col = hcl.colors(100, "Blues", rev = TRUE))
  
  return(stats)
}

# Process each chunk file sequentially
results <- lapply(chunk_files, process_chunk)

# Combine all results into a single dataframe
count_obs_df <- do.call(rbind, results)

# View the results
print(count_obs_df)
