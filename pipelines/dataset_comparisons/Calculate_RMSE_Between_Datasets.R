library(terra)

# Load your raster stacks
amsr <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-Unified_12km/stack/AMSR-Unified_12km_Harmonization.nc")
nsidc <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/25km_Sea-Ice-Index/stack/NSIDC-Sea-Ice-Index_25km_Continental_Harmonization.nc")

# Apply the threshold to remove values from 0 to 15
threshold <- function(x) {
  x[x <= 15] <- NA
  return(x)
}

amsr <- app(amsr, threshold)
nsidc <- app(nsidc, threshold)


# Compute squared differences for each cell across all layers
squared_diffs <- (amsr - nsidc)^2

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
plot(cell_rmse, main="Cell-level RMSE", colNA="blue")