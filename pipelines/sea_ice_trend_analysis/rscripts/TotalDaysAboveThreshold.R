library(terra)
library(lubridate)

# raster1 <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/raster-stack/Chunk_2014-01-21.nc")

# # Assuming raster1 is already loaded and set
raster_dates <- as.Date(time(raster1))
July_2014_layers <- raster1[[raster_dates %within% interval("2014-01-01", "2014-12-31")]]


# Define a threshold
threshold <- 50

# Create a binary stack where each cell is 1 if above threshold, otherwise 0
binary_stack <- July_2014_layers > threshold


# Correctly sum the binary stack to count days above threshold per cell
# Instead of sum(), we use the app() function to apply a sum across layers
count_days_above_threshold <- app(binary_stack, sum)

# Plot the result
plot(count_days_above_threshold, main="Number of Days Above Threshold")



