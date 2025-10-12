# Load the necessary library
library(dplyr)

# Load the dataset
file_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/AllCounts_V_4_1.csv"
data <- read.csv(file_path)

# Filter for gentoo penguins only in CAMLR region 48.1
gentoo_data <- data %>%
  filter(tolower(common_name) == "gentoo penguin" & cammlr_region == 48.1)

# Export the filtered data to a new CSV file
write.csv(gentoo_data, "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/GentooCounts_48_1.csv", row.names = FALSE)

# Print a message indicating the file has been saved
cat("Filtered data has been saved to GentooCounts_48_1.csv")
