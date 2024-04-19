library(terra)

# Function to extract the date from a filename
extractDate <- function(filename) {
  base_name <- basename(filename)
  date_string <- sub("S_(\\d{8})_.*", "\\1", base_name)
  date_converted <- as.Date(date_string, format="%Y%m%d")
  if (is.na(date_converted)) {
    warning("Date conversion failed for: ", base_name)
  }
  return(date_converted)
}

# Define the directory with raster files
raster_dir <- "C:/Users/wethi/Desktop/dwnld/NSIDC-Sea-Ice-Index/tif"
raster_files <- list.files(path = raster_dir, pattern = "\\.tif$", full.names = TRUE)
raster_files <- raster_files[order(sapply(raster_files, extractDate))]

# Calculate the size of each chunk
num_files <- length(raster_files)
chunk_size <- ceiling(num_files / 4)  # Adjust the denominator to change chunking strategy

# Function to process each chunk
process_chunk <- function(files, chunk_index) {
  print(paste("Processing chunk", chunk_index, "with", length(files), "files..."))
  raster_stack <- rast(files)
  names(raster_stack) <- sapply(files, function(x) format(extractDate(x), "%Y-%m-%d"))
  dates <- as.Date(sapply(files, extractDate))
  time(raster_stack) <- dates
  values_to_exclude <- c(0, 2510, 2530, 2540, 2550)
  raster_stack[raster_stack %in% values_to_exclude] <- NA
  raster_stack[raster_stack < 150] <- NA
  raster_stack <- raster_stack / 10
  
  # Define output paths for NetCDF
  nc_output_path <- sprintf("D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/raster-stack/Chunk_%s.nc", format(dates[1], "%Y-%m-%d"))
  
  # Write to NetCDF
  writeCDF(raster_stack, filename = nc_output_path, varname="sea_ice_concentration", longname="Sea Ice Concentration", unit="percentage", overwrite=TRUE, compression=9, missval=-9999)
  print(paste("Chunk", chunk_index, "processed and saved to", nc_output_path))
  return(nc_output_path)
}

# Loop through each chunk
for (i in seq(1, num_files, by = chunk_size)) {
  chunk_index <- (i - 1) / chunk_size + 1
  chunk_files <- raster_files[i:min(i + chunk_size - 1, num_files)]
  output_path <- process_chunk(chunk_files, chunk_index)
  
  # Clear memory
  rm(raster_stack)
  gc()  # Force garbage collection
  print(paste("Memory cleared after processing chunk", chunk_index))
}
v
# Optionally, load and verify one of the outputs
if (file.exists(output_path)) {
  sample_output <- rast(output_path)
  print(sample_output)
  print(time(sample_output))
} else {
  print(paste("No output at path:", output_path))
}
