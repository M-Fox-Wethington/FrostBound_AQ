library(terra)

# # Load your raster stacks
amsr_monthly <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/12km_AMSR-Unified/stack/substack/AMSR-Unified_harmonization_winter_monthly_mean_extent.nc")
nsidc_monthly <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/25km_Sea-Ice-Index/stack/substack/25km_Sea-Ice-Index_harmonization_winter_monthly_mean_extent.nc")


# Apply the threshold to remove values from 0 to 15
threshold <- function(x) {
  x[x <= 15] <- NA
  return(x)
}

amsr_monthly <- app(amsr_monthly, threshold)
nsidc_monthly <- app(nsidc_monthly, threshold)

# Compute squared differences for each cell across all layers
squared_diffs <- (amsr_monthly - nsidc_monthly)^2

# Function to calculate RMSE across layers
calc_rmse <- function(x) {
  n <- sum(!is.na(x))
  if (n > 0) {
    mse <- sum(x, na.rm = TRUE) / n
    return(sqrt(mse))
  } else {
    return(NA)
  }
}

# Apply the RMSE function across all cells
cell_rmse <- app(squared_diffs, calc_rmse)

# Plot the cell-level RMSE raster
plot(cell_rmse, main="Cell-level RMSE - Monthly Averages", colNA="blue")

# Calculate and print the average and median RMSE
average_rmse <- mean(values(cell_rmse), na.rm = TRUE)
median_rmse <- median(values(cell_rmse), na.rm = TRUE)

# Print the results
print(paste("Average RMSE: ", average_rmse))
print(paste("Median RMSE: ", median_rmse))
