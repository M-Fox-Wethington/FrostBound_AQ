
# Load necessary libraries
library(ggplot2)
library(hexbin)
library(viridis)

# Function to calculate RMSE (Root Mean Squared Error)
calc_rmse <- function(actual, predicted) {
  return(sqrt(mean((actual - predicted)^2)))
}

# Extract the hexbin cell IDs and counts
hb_cID <- hb@cell  # This contains the hexbin IDs for each data point

# Initialize a vector to store RMSE values for each hexbin
rmses <- numeric(length(hb@count))  # Use hb@count, which stores the counts for each hexbin

# Loop through each hexbin and calculate RMSE
for (i in seq_along(hb@cell)) {
  # Get the indices of the points that fall in the i-th hexbin
  bin_indices <- which(hb_cID == i)
  
  # Get the AMSR and NSIDC values for these points
  if (length(bin_indices) > 0) {
    amsr_values <- data$AMSR[bin_indices]
    nsidc_values <- data$NSIDC[bin_indices]
    
    # Calculate RMSE for this bin
    rmses[i] <- calc_rmse(amsr_values, nsidc_values)
  } else {
    rmses[i] <- NA  # Assign NA if there are no values in the bin
  }
}

# Now, create the data frame with RMSE and hexbin coordinates
hb_data <- data.frame(hcell2xy(hb), count = hb@count, RMSE = rmses)

# Remove rows with NA RMSE values
hb_data <- hb_data[!is.na(hb_data$RMSE), ]

# Plot the hexbin plot with RMSE as the color
ggplot(hb_data, aes(x = x, y = y, fill = RMSE)) +
  geom_hex(stat = "identity") +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "RMSE") +
  labs(
    title = "Root Mean Squared Error (RMSE) - AMSR vs. NSIDC",
    x = "AMSR Unified 12.5 km (Actual, 0-1)",
    y = "NSIDC 12.5 km Harmonized (Predicted, 0-1)"
  ) +
  theme_minimal()
