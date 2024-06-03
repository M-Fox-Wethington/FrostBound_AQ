library(tidyverse)


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

# create site x season template which is used throughout analysis
w_template <- SiteList %>%
  dplyr::select(site_id, site) %>%
  # expand each site by the number of seasons
  uncount(n_seasons) %>%
  # create relative season index for each site
  mutate(season_relative = rep(1:n_seasons, n_sites)) %>%
  # create season var from relative season index
  mutate(season = season_relative + min_season - 1) %>%
  arrange(season_relative, site)

w_df <- rbind(
  # keep all presence/absence data and assign observation type as 2 (observed)
  mapppdr::penguin_obs %>%
    dplyr::filter(species_id == species & season >= min_season & season <= max_season) %>%
    dplyr::select(site_id, season, presence) %>%
    mutate(known_w = 1),
  # append presence/absence assumption data which is not part of mapppd
  # and assign observation type of 1 (assumed)
  data.frame(read_csv(file = "D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/gentoo_presence_absence_assumptions.csv")) %>%
    mutate(known_w = 0)) %>%
  # determine for each site x season if breeding is observed or assumed
  group_by(site_id, season) %>%
  summarise(w = base::max(presence), known_w = base::max(known_w)) %>%
  ungroup() %>%
  # join with w_template to fill in missing site x seasons with no presence/absence data
  right_join(w_template, by = c("site_id", "season")) %>%
  # assign observation type as 0 (imputed)
  mutate(known_w = replace(known_w, is.na(known_w), 0)) %>%
  arrange(site_id, season) %>%
  # impute missing presence/absence data using the following assumptions
  # ASSSUMPTION: fill in NA between (1,1) with 1
  # ASSSUMPTION: fill in NA between (0,1) with 0
  # ASSSUMPTION: fill in NA between (1,0) with 1
  # ASSSUMPTION: fill in NA between (.,1) and (1,.) with 1
  # ASSSUMPTION: fill in NA between (.,0) and (0,.) with 0
  dplyr::group_by(site_id) %>%
  tidyr::fill(w, .direction = "downup") %>%
  dplyr::ungroup() %>%
  # create second site_id var for plotting sites alphabetically in ggplot
  mutate(site_id_rev = factor(site_id, levels = rev(sort(unique(site_id))))) %>%
  dplyr::select(site_id, site_id_rev, season, site, season_relative, w, known_w)

# convert w to matrix to be used in model
w <- w_df %>%
  dplyr::select(site, season_relative, w) %>%
  # create matrix where rows are sites and columns are seasons
  pivot_wider(names_from = season_relative, values_from = w, names_sort = TRUE) %>%
  dplyr::select(-site) %>%
  as.matrix()


abundance <- mapppdr::penguin_obs %>%
  # keep all counts between min and max season
  dplyr::filter(count > 0 & species_id == species & season >= min_season & season <= max_season) %>%
  # join to get site index and initial season
  right_join(SiteList, by = "site_id") %>%
  # create relative season index 
  mutate(season_relative = season - min_season + 1) %>%
  # ASSUMPTION: increase accuracy category of all adult counts by + 3 with a max error of 5
  rowwise() %>%
  mutate(accuracy = replace(accuracy, type == "adults", base::min((accuracy[type == "adults"] + 3), 5))) %>%
  ungroup() %>%  
  mutate(type = replace(type, type == "adults", "nests")) %>%
  # ASSUMPTION: keep maximum nest and chick count reported each season for a site
  group_by(site_id, season, season_relative, type) %>%
  arrange(desc(count), accuracy) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  # ASSUMPTION: convert accuracy to the following errors/precisions
  mutate(sigma = case_when(
    accuracy == 1 ~ 0.02490061, 
    accuracy == 2 ~ 0.04955838,
    accuracy == 3 ~ 0.1201131, 
    accuracy == 4 ~ 0.2212992, 
    accuracy == 5 ~ 0.4472728)) %>%
  mutate(precision = case_when(
    accuracy == 1 ~ 1/0.02490061^2, 
    accuracy == 2 ~ 1/0.04955838^2,
    accuracy == 3 ~ 1/0.1201131^2, 
    accuracy == 4 ~ 1/0.2212992^2, 
    accuracy == 5 ~ 1/0.4472728^2)) %>%  
  dplyr::select(site_id, site, season, season_relative, initial_season, type, 
                count, presence, accuracy, sigma, precision) %>%
  arrange(site, season_relative, type, -count, accuracy, sigma, precision)  

abundance_initial <- abundance %>%
  # keep first observed count for each site's time series
  dplyr::filter(initial_season == season_relative) %>%
  # ASSUMPTION: if no nest count is available in the initial season and a chick count is then
  # assume chick count is 1:1 nest count
  group_by(site_id, season, site, season_relative) %>%
  arrange(desc(type)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(site_id, season, site, season_relative, count, sigma, precision)

abundance_nests <- abundance %>%
  # keep all nest counts after the initial season
  dplyr::filter(initial_season != season_relative & type == "nests") %>%
  dplyr::select(site_id, season, site, season_relative, count, sigma, precision)

abundance_chicks <- rbind(
  # keep all chick counts after the initial season
  abundance %>%
    dplyr::filter(initial_season != season_relative & type == "chicks") %>%
    dplyr::select(site_id, season, site, season_relative, count, sigma, precision),
  # append chick counts from initial season that were not converted to nest counts
  # meaning there was both a chick and nest count in the initial season
  abundance %>%
    dplyr::filter(initial_season == season_relative) %>%
    group_by(site_id, season, site, season_relative) %>%
    arrange(desc(type)) %>%
    slice(2) %>%
    ungroup() %>%
    dplyr::select(site_id, season, site, season_relative, count, sigma, precision))

# moment match alpha shape and rate parameters for breeding productivity
mu <- .5
sigma <- .25
a <- (mu^2 - mu^3 - mu * sigma^3) / sigma^2
b <- (mu - 2* mu^2 + mu^3 - sigma^2 + mu * sigma^3) / sigma^2

# create the data list for the JAGS model
model_data <- list(
  nests = nrow(abundance_nests),
  y_n = log(abundance_nests$count), 
  precision_n = abundance_nests$precision,
  site_n = abundance_nests$site,
  season_n = abundance_nests$season_relative,
  chicks = nrow(abundance_chicks),
  y_c = log(abundance_chicks$count), 
  precision_c = abundance_chicks$precision,
  site_c = abundance_chicks$site,
  season_c = abundance_chicks$season_relative,
  y_i = log(abundance_initial$count),
  precision_i = abundance_initial$precision,
  n_sites = n_sites,
  n_seasons = n_seasons,
  s = as.vector(SiteList$initial_season),
  w = w,
  a = a,
  b = b)
