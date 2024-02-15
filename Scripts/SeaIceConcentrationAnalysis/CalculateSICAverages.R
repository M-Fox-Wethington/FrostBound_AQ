# Define a function to calculate sea ice concentration (SIC) averages for a specified region and time period
# using satellite images.
# Args:
#   ssmi_tif_dir: Directory containing TIFF files of SSMI data
#   colony_buffer_shp_dir: Directory of the shapefile representing the buffer around a colony
#   region: The region of interest as a string
#   time_period: The specific time period (e.g., Annual, Q1, Q2, Q3, Q4) for analysis
#   out_dir: Output directory for the resulting CSV file

ssmi_SeaIceAverages <- function(ssmi_tif_dir, colony_buffer_shp_dir, region, time_period, out_dir){
  
  # Load necessary libraries for data manipulation, spatial analysis, and plotting
  require(terra) # For spatial data analysis and raster manipulation
  require(sf) # For handling vector data
  require(data.table) # For efficient data manipulation
  require(tidyverse) # For data manipulation and visualization
  require(tibble) # For creating tibbles, a modern reimagining of data frames
  require(stringr) # For string manipulation
  require(ggplot2) # For data visualization
  
  # Configure R to use fixed scientific notation for readability
  options("scipen"= 100, "digits"=4)
  
  # Generate a list of all TIFF files within the specified directory
  rast_list <- list.files(ssmi_tif_dir, 
                          pattern='.tif$', 
                          all.files= TRUE, 
                          full.names=TRUE)
  
  # Read the shapefile representing the buffer zone around a colony
  v <- vect(colony_buffer_shp_dir)
  
  # Create a raster stack from the list of TIFF files
  r <- rast(rast_list)
  
  # Clip the raster to the extent of the buffer zone plus a small margin (0.01 degrees)
  crop <- terra::crop(r, ext(v) + .01)
  
  # Mask the clipped raster with the buffer zone to isolate the area of interest
  mask <- mask(crop, v)
  
  # Convert the pixel values to represent SIC percentages (original scale is 0-1000)
  data.SIC <- (mask / 10) 
  
  # Calculate the global mean SIC, removing NA values
  data.SIC.frame <-  as.data.frame(terra::global(data.SIC, 
                                                 fun = mean, 
                                                 na.rm=TRUE))
  
  # Convert row names to a separate column for further manipulation
  data.SIC.frame <- tibble::rownames_to_column(data.SIC.frame, "rID")
  
  # Duplicate the rID column to extract date information
  data.SIC.frame$Date <- data.SIC.frame$rID
  
  # Extract the date substring and format it for readability
  data.SIC.frame <- data.SIC.frame %>%
    mutate(Date = substr(Date, 3, 8)) %>% 
    rename(Mean_SIC = mean)
  
  # Reformat the date to include a dash for better readability
  data.SIC.frame$Date <- sub("(.{4})(.*)", "\\1-\\2", data.SIC.frame$Date)
  
  # Split the formatted date into separate Year and Month columns
  data.SIC.frame <- data.SIC.frame %>% separate(Date, c('Year', 'Month'))
  
  # Calculate average SIC based on the specified time period (Annual, Q1, Q2, Q3, Q4)
  # and group data accordingly for analysis
  if(time_period == "Annual")
  { 
    season <- "Annual" 
    data.SIC.frame <- data.SIC.frame %>%
      group_by(Year) %>%
      summarize(GlobalMean = mean(Mean_SIC, na.rm=TRUE)) 
  }
  
  # Repeat the process for each quarter, filtering by the respective months
  # This approach isolates the data relevant to each season before calculating the average SIC
  if(time_period == "Q1") { 
    season <- "Q1" 
    # Filter for Q1 months (January, February, March)
    data.SIC.frame <- data.SIC.frame %>%
      filter(Month %in% c("01", "02", "03")) %>%
      group_by(Year) %>%
      summarize(GrandMean = mean(Mean_SIC, na.rm=TRUE)) 
  }
  
  # Similar filtering and summarization for Q2, Q3, and Q4
  if(time_period == "Q2"){  
    season <- "Q2" 
    data.SIC.frame <- data.SIC.frame %>%
      filter(Month %in% c("04", "05", "06")) %>%
      group_by(Year) %>%
      summarize(GrandMean = mean(Mean_SIC, na.rm=TRUE))
  }
  
  if(time_period == "Q3"){  
    season <- "Q3" 
    data.SIC.frame <- data.SIC.frame %>%
      filter(Month %in% c("07", "08", "09")) %>%
      group_by(Year) %>%
      summarize(GrandMean = mean(Mean_SIC, na.rm=TRUE))
  }
  
  if(time_period == "Q4"){  
    season <- "Q4" 
    data.SIC.frame <- data.SIC.frame %>%
      filter(Month %in% c("10", "11", "12")) %>%
      group_by(Year) %>%
      summarize(GrandMean = mean(Mean_SIC, na.rm=TRUE))
  }
  
  # Optionally, view the head of the dataframe to inspect the first few rows
  # This line can be commented out or removed based on preference
  head(data.SIC.frame)
  
  # Construct the output file name dynamically based on the region and season
  output_file_name <- paste0(out_dir, "/", region, "_", season, "_Mean_SIC_1979-2021.csv")
  
  # Export the data frame to a CSV file in the specified output directory
  write.csv(data.SIC.frame, output_file_name, row.names = FALSE)
  
  # Print a message to indicate successful completion and the location of the output file
  message("Data exported successfully to ", output_file_name)
}
