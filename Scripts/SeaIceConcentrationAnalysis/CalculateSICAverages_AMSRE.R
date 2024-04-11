library(terra)    # For spatial data analysis and raster manipulation
library(dplyr)    # For data manipulation and analysis
library(tibble)   # For creating and manipulating tibbles
library(stringr)  # For string manipulation
library(lubridate) # For date manipulation

calculateSICaveragesForRegion <- function(raster_stack, region_shp_path, out_dir) {
  # Read the shapefile representing the region of interest
  region <- vect(region_shp_path)
  
  # The target CRS will be the CRS of the raster_stack
  target_crs <- crs(raster_stack)
  
  # Check if the shapefile's CRS matches the raster stack's CRS; if not, reproject
  if (crs(region) != target_crs) {
    region <- project(region, target_crs)
  }
  
  masked_stack <- mask(raster_stack, region)
  
  # Convert the pixel values to represent SIC percentages (original scale is 0-1000)
  masked_stack[masked_stack > 100] <- NA
  
  # Initialize a list to store SIC data frames
  sic_list <- list()
  
  # Extract layer names, assumed to be dates in YYYYMMDD format
  layer_dates <- names(masked_stack)
  
  print(layer_dates)
  
  # Loop through each layer to calculate SIC
  for (date in layer_dates) {
    parsed_date <- ymd(date)
    year <- year(parsed_date)
    month <- month(parsed_date)
    week <- week(parsed_date)
    
    mean_sic <- global(masked_stack[[date]], fun = 'mean', na.rm = TRUE)
    mean_sic <- as.numeric(mean_sic) # Ensure mean_sic is numeric
    
    print(paste0("Mean SIC for Current layer: ", mean_sic))
    
    sic_df <- tibble(Date = parsed_date, Year = year, Month = month, Week = week, Mean_SIC = mean_sic)
    sic_list[[date]] <- sic_df
  }
  
  
  all_sic_data <- bind_rows(sic_list)
  
  # Now you can check for any non-numeric values in all_sic_data$Mean_SIC
  if(any(!is.numeric(all_sic_data$Mean_SIC))) {
    stop("Non-numeric values found in Mean_SIC.")
  }
  
  # Combine all daily data frames into one
  all_sic_data <- bind_rows(sic_list)
  
  # Calculate weekly, monthly, and annual averages
  weekly_sic <- all_sic_data %>%
    group_by(Year, Week) %>%
    summarize(WeeklyMean_SIC = mean(Mean_SIC, na.rm = TRUE), .groups = 'drop')
  
  monthly_sic <- all_sic_data %>%
    group_by(Year, Month) %>%
    summarize(MonthlyMean_SIC = mean(Mean_SIC, na.rm = TRUE), .groups = 'drop')
  
  annual_sic <- all_sic_data %>%
    group_by(Year) %>%
    summarize(AnnualMean_SIC = mean(Mean_SIC, na.rm = TRUE), .groups = 'drop')
  
  # Export the data to CSV files
  write.csv(weekly_sic, file.path(out_dir, "SIC_weekly_averages.csv"), row.names = FALSE)
  write.csv(monthly_sic, file.path(out_dir, "SIC_monthly_averages.csv"), row.names = FALSE)
  write.csv(annual_sic, file.path(out_dir, "SIC_annual_averages.csv"), row.names = FALSE)
  
  # Print success messages
  message("Weekly, monthly, and annual SIC averages successfully calculated and saved to ", out_dir)
}



# Example usage of the function:
calculateSICaveragesForRegion(
  raster_stack = loaded_raster_stack,
  region_shp_path = "D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/CCAMLR_MPA_Planning_Domains_WAP.shp",
  out_dir = "D:/Manuscripts_localData/FrostBound_AQ/Results"
)