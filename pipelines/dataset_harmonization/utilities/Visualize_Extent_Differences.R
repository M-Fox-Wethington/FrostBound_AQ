# Load necessary libraries
library(terra)
library(ggplot2)
library(patchwork)

# Load the raster black mask (we will color all masked areas black)
black_mask_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/Land_Mask/Black_Mask.tif"
black_mask <- rast(black_mask_path)

# Function to calculate the pixel-level difference (retain direction)
calculate_pixel_diff_raster <- function(layer1, layer2) {
  diff <- layer1 - layer2  # Calculate the difference (retain direction)
  return(diff)  # Return the difference as a SpatRaster object
}

# Select dates for the layers to be compared (adjust according to your needs)
date1 <- as.Date(amsr_dates[1])
date2 <- as.Date(amsr_dates[2])
date3 <- as.Date(amsr_dates[3])

# Calculate pixel-level differences for AMSR vs. NSIDC (original)
pixel_diff_amsr_nsidc_1 <- calculate_pixel_diff_raster(amsr_stack[[1]], nsidc_stack[[1]])
pixel_diff_amsr_nsidc_2 <- calculate_pixel_diff_raster(amsr_stack[[2]], nsidc_stack[[2]])
pixel_diff_amsr_nsidc_3 <- calculate_pixel_diff_raster(amsr_stack[[3]], nsidc_stack[[3]])

# Calculate pixel-level differences for AMSR vs. Harmonized
pixel_diff_amsr_harmonized_1 <- calculate_pixel_diff_raster(amsr_stack[[1]], nsidc_harmonized_stack[[1]])
pixel_diff_amsr_harmonized_2 <- calculate_pixel_diff_raster(amsr_stack[[2]], nsidc_harmonized_stack[[2]])
pixel_diff_amsr_harmonized_3 <- calculate_pixel_diff_raster(amsr_stack[[3]], nsidc_harmonized_stack[[3]])

# Convert the raster data to data frames for plotting
df_amsr_nsidc_1 <- as.data.frame(pixel_diff_amsr_nsidc_1, xy = TRUE)
df_amsr_nsidc_2 <- as.data.frame(pixel_diff_amsr_nsidc_2, xy = TRUE)
df_amsr_nsidc_3 <- as.data.frame(pixel_diff_amsr_nsidc_3, xy = TRUE)

df_amsr_harmonized_1 <- as.data.frame(pixel_diff_amsr_harmonized_1, xy = TRUE)
df_amsr_harmonized_2 <- as.data.frame(pixel_diff_amsr_harmonized_2, xy = TRUE)
df_amsr_harmonized_3 <- as.data.frame(pixel_diff_amsr_harmonized_3, xy = TRUE)

# Convert black mask to a data frame
df_black_mask <- as.data.frame(black_mask, xy = TRUE)

# Define custom diverging color scale
custom_colors <- scale_fill_gradient2(
  low = "steelblue",    # Negative differences (blue)
  mid = "white",        # Zero difference (white)
  high = "orange",      # Positive differences (orange)
  na.value = "black",   # NA values as black for background mask
  limits = c(-0.5, 0.5),  # Adjust the range as needed
  name = "Difference"
)

# Plot AMSR vs. NSIDC (original) pixel differences with black mask as background
plot_amsr_nsidc_1 <- ggplot() +
  geom_raster(data = df_black_mask, aes(x = x, y = y), fill = "black", na.rm = TRUE) +  # Plot black mask first (background)
  geom_raster(data = df_amsr_nsidc_1, aes(x = x, y = y, fill = df_amsr_nsidc_1[, 3])) +
  custom_colors +
  labs(title = paste("AMSR vs NSIDC (Original) - Date:", date1)) +
  theme_minimal()

plot_amsr_nsidc_2 <- ggplot() +
  geom_raster(data = df_black_mask, aes(x = x, y = y), fill = "black", na.rm = TRUE) +  # Plot black mask first (background)
  geom_raster(data = df_amsr_nsidc_2, aes(x = x, y = y, fill = df_amsr_nsidc_2[, 3])) +
  custom_colors +
  labs(title = paste("AMSR vs NSIDC (Original) - Date:", date2)) +
  theme_minimal()

plot_amsr_nsidc_3 <- ggplot() +
  geom_raster(data = df_black_mask, aes(x = x, y = y), fill = "black", na.rm = TRUE) +  # Plot black mask first (background)
  geom_raster(data = df_amsr_nsidc_3, aes(x = x, y = y, fill = df_amsr_nsidc_3[, 3])) +
  custom_colors +
  labs(title = paste("AMSR vs NSIDC (Original) - Date:", date3)) +
  theme_minimal()

# Plot AMSR vs. Harmonized pixel differences with black mask as background
plot_amsr_harmonized_1 <- ggplot() +
  geom_raster(data = df_black_mask, aes(x = x, y = y), fill = "black", na.rm = TRUE) +  # Plot black mask first (background)
  geom_raster(data = df_amsr_harmonized_1, aes(x = x, y = y, fill = df_amsr_harmonized_1[, 3])) +
  custom_colors +
  labs(title = paste("AMSR vs Harmonized - Date:", date1)) +
  theme_minimal()

plot_amsr_harmonized_2 <- ggplot() +
  geom_raster(data = df_black_mask, aes(x = x, y = y), fill = "black", na.rm = TRUE) +  # Plot black mask first (background)
  geom_raster(data = df_amsr_harmonized_2, aes(x = x, y = y, fill = df_amsr_harmonized_2[, 3])) +
  custom_colors +
  labs(title = paste("AMSR vs Harmonized - Date:", date2)) +
  theme_minimal()

plot_amsr_harmonized_3 <- ggplot() +
  geom_raster(data = df_black_mask, aes(x = x, y = y), fill = "black", na.rm = TRUE) +  # Plot black mask first (background)
  geom_raster(data = df_amsr_harmonized_3, aes(x = x, y = y, fill = df_amsr_harmonized_3[, 3])) +
  custom_colors +
  labs(title = paste("AMSR vs Harmonized - Date:", date3)) +
  theme_minimal()

# Combine the plots side by side for AMSR vs NSIDC and AMSR vs Harmonized comparisons
combined_plot <- (plot_amsr_nsidc_1 + plot_amsr_harmonized_1) / 
  (plot_amsr_nsidc_2 + plot_amsr_harmonized_2) / 
  (plot_amsr_nsidc_3 + plot_amsr_harmonized_3) +
  plot_layout(guides = 'collect')

# Print the combined plot
print(combined_plot)
