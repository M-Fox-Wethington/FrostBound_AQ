library(terra)
library(sf)

# Directory and files setup
input_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/NASA_Bootstrap_sic/staged"
output_dir <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/NASA_Bootstrap_sic/processed"

# Initialize lists to store filenames that weren't successfully processed
unprocessed_files <- list()
files_without_icecon <- list()

# Define a function to process a single NetCDF file
process_netcdf_file <- function(nc_file, input_dir, output_dir) {
  tryCatch({
    r <- rast(nc_file)
    
    # Print all variable names for diagnostic purposes
    print(names(r))
    
    # Find any variables that contain 'ICECON' in their name
    available_vars <- grep("ICECON", names(r), value = TRUE)
    
    if (length(available_vars) > 0) {
      sic_raster <- rast(nc_file, subds=available_vars[1])
      
      values(sic_raster)[values(sic_raster) == 1100] <- NA  # Missing data
      values(sic_raster)[values(sic_raster) == 1200] <- NA  # Land or ice-shelf
      values(sic_raster)[values(sic_raster) == 0] <- NA  # Land or ice-shelf
      
      # Use the filename from the input for the output
      output_filename <- gsub(input_dir, output_dir, nc_file)
      output_filename <- gsub(".nc", ".tif", output_filename)
      
      # Save the raster data as a GeoTIFF
      writeRaster(sic_raster, filename=output_filename, overwrite=TRUE)
      print(paste0("Raster saved as: ", output_filename))
    } else {
      files_without_icecon <<- append(files_without_icecon, nc_file)
      print(paste0("No ICECON data in ", nc_file))
    }
  }, error = function(e) {
    unprocessed_files <<- append(unprocessed_files, nc_file)
    print(paste0("Error processing ", nc_file, ": ", e$message))
  })
}




# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Process each NetCDF file
netcdf_files <- list.files(input_dir, pattern = "\\.nc$", full.names = TRUE)
lapply(netcdf_files, function(file) process_netcdf_file(file, input_dir, output_dir))
# After the lapply call
unprocessed_files_df <- data.frame(files = unlist(unprocessed_files), stringsAsFactors = FALSE)
files_without_icecon_df <- data.frame(files = unlist(files_without_icecon), stringsAsFactors = FALSE)

# Write to CSV
write.csv(unprocessed_files_df, file=paste0(output_dir, "unprocessed_files.csv"), row.names = FALSE)
write.csv(files_without_icecon_df, file=paste0(output_dir, "files_without_icecon.csv"), row.names = FALSE)


