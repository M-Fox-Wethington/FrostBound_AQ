library(tidyverse)
library(coda)
library(dplyr)
library(tidyr)

required_packages <- c("rlang", "fastmap", "digest", "fs", "cachem", "vctrs", "stringi", "glue", "cli", "utf8", "fansi", "dplyr")

for(pkg in required_packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
  }
}

min_season <- 1970
max_season <- 2023
species <- "GEPE"

# assign the total number of seasons as n_seasons
(n_seasons <- (max_season - min_season) + 1)

# Define SiteList using mapppdr data (assumed to be loaded previously)
SiteList <- mapppdr::penguin_obs %>%
  # keep all sites that have at least 1 count between min and max season
  dplyr::filter(count > 0 & species_id == species & season >= min_season & season <= max_season) %>%
  # create relative season index 
  mutate(season_relative = season - min_season + 1) %>%
  # determine first season a count is observed for each site
  group_by(site_id) %>%
  summarise(initial_season = min(season_relative)) %>%
  ungroup() %>%
  # join to get other site specific covariates for visualization purposes
  left_join(mapppdr::sites, by = "site_id") %>%
  # create site index for model and visualization
  mutate(site = as.numeric(as.factor(site_id))) %>%
  dplyr::select(site_id, site_name, ccamlr_id, site, initial_season, latitude, longitude)

(n_sites <- nrow(SiteList))

# Load the MCMC samples
load("D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/model_data_rinits_output.rda")

# Convert the MCMC samples to a matrix
mcmc_matrix <- as.matrix(model_data_rinits_output)

# Extract the growth rate estimates
growth_rate_samples <- mcmc_matrix[, grep("^l_a\\[", colnames(mcmc_matrix))]

# Check the structure of the extracted growth rate samples
str(growth_rate_samples)

# Updated extract_indices function
extract_indices <- function(colname) {
  indices <- gsub("[^0-9,]", "", colname)
  as.integer(unlist(strsplit(indices, ",")))
}

indices <- lapply(colnames(growth_rate_samples), extract_indices)
sites <- sapply(indices, `[`, 1)
seasons <- sapply(indices, `[`, 2)

# Calculate the mean and credible intervals for each site and season
growth_rate_summary <- apply(growth_rate_samples, 2, function(x) {
  c(mean = mean(x), 
    lower_95 = quantile(x, 0.025), 
    upper_95 = quantile(x, 0.975))
})

# Convert the summary matrix to a dataframe
growth_rate_df <- as.data.frame(t(growth_rate_summary))
growth_rate_df$site <- sites
growth_rate_df$season <- seasons

# Add the year column based on the season and min_season
growth_rate_df <- growth_rate_df %>%
  mutate(year = season + min_season - 1)

# Arrange the dataframe by site and season
growth_rate_df <- growth_rate_df %>%
  arrange(site, season)

# View the dataframe
print(growth_rate_df)

# Perform the left join to include additional site-specific covariates
growth_rate_df_with_site_info <- growth_rate_df %>%
  left_join(SiteList %>% select(site, site_id, site_name, ccamlr_id, latitude, longitude), by = "site")

# Print the resulting data frame to verify the join
print(growth_rate_df_with_site_info)
str(growth_rate_df_with_site_info)
