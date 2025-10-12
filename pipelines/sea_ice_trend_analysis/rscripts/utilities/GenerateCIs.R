# Load necessary libraries
library(dplyr)

# Load your data from CSV
data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/trend-analysis/Trends/significant_results.csv")

# Calculate the 95% CI for the Intercept and Time estimates
data_with_cis <- data %>%
  mutate(
    Intercept_CI_Lower = Intercept_Estimate - 1.96 * Intercept_StdError,
    Intercept_CI_Upper = Intercept_Estimate + 1.96 * Intercept_StdError,
    Time_CI_Lower = Time_Estimate - 1.96 * Time_StdError,
    Time_CI_Upper = Time_Estimate + 1.96 * Time_StdError
  )

# Save the updated data with CIs to a new CSV file
output_path <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/trend-analysis/Trends/significant_results_with_CIs.csv"
write.csv(data_with_cis, output_path, row.names = FALSE)

# View the first few rows of the updated data
head(data_with_cis)




# Load necessary libraries
library(dplyr)

# Load your data from the specified CSV file
data_2010_2024 <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/trend-analysis/Trends/2010-2024/significant_results.csv")

# Calculate the 95% CI for the Intercept and Time estimates
data_with_cis_2010_2024 <- data_2010_2024 %>%
  mutate(
    Intercept_CI_Lower = Intercept_Estimate - 1.96 * Intercept_StdError,
    Intercept_CI_Upper = Intercept_Estimate + 1.96 * Intercept_StdError,
    Time_CI_Lower = Time_Estimate - 1.96 * Time_StdError,
    Time_CI_Upper = Time_Estimate + 1.96 * Time_StdError
  )

# Save the updated data with CIs to a new CSV file
output_path_2010_2024 <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/dataset-harmonization/trend-analysis/Trends/2010-2024/significant_results_with_CIs.csv"
write.csv(data_with_cis_2010_2024, output_path_2010_2024, row.names = FALSE)

# View the first few rows of the updated data
head(data_with_cis_2010_2024)
