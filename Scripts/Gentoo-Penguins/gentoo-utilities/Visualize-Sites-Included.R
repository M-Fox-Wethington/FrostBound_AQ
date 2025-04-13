# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load Gentoo penguin data
penguin_data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Datasets/mapppd/gentoo_presence_absence_assumptions.csv")
penguin_abundance_data <- read.csv("D:/Manuscripts_localData/FrostBound_AQ/Datasets/gentoo-abundance-model/inputs/modeled_gentoo_parameters.csv")


# Adjust the season column in the penguin abundance data
penguin_abundance_data <- penguin_abundance_data %>%
  mutate(year = 1970 + season - 1)

# Filter penguin data to keep only consecutively censused sites and clean duplicates
consecutively_censused <- penguin_data %>%
  arrange(site_id, season) %>%
  group_by(site_id) %>%
  mutate(next_season = lead(season),
         next_presence = lead(presence)) %>%
  filter(presence == 1 & next_presence == 1 & season == next_season - 1) %>%
  ungroup() %>%
  select(site_id, season, presence, next_season)

# Merge with penguin abundance data
penguin_abundance_filtered <- penguin_abundance_data %>%
  left_join(consecutively_censused, by = c("site_id", "year" = "season")) %>%
  filter(!is.na(next_season))

# Ensure first year in consecutive pairs is also flagged
penguin_abundance_filtered <- penguin_abundance_filtered %>%
  group_by(site_id) %>%
  mutate(consecutive_used = ifelse(presence == 1 & lag(presence, default = 0) == 1, "Used", "Not Used")) %>%
  ungroup()

penguin_abundance_filtered <- penguin_abundance_filtered %>%
  group_by(site_id) %>%
  mutate(consecutive_used = ifelse(lead(consecutive_used, default = "Not Used") == "Used", "Used", consecutive_used)) %>%
  ungroup()

# Create a complete grid of site_id and year to fill in missing values
site_year_grid <- expand.grid(
  site_id = unique(penguin_abundance_filtered$site_id),
  year = seq(min(penguin_abundance_filtered$year), max(penguin_abundance_filtered$year))
)

# Merge with the filtered data to identify missing years
penguin_abundance_complete <- site_year_grid %>%
  left_join(penguin_abundance_filtered, by = c("site_id", "year")) %>%
  mutate(presence_indicator = ifelse(is.na(presence), "Absent", "Observed"),
         consecutive_used = ifelse(is.na(consecutive_used), "Not Used", consecutive_used))

# Filter to keep only the rows that were observed and used in the model
penguin_abundance_used <- penguin_abundance_complete %>%
  filter(consecutive_used == "Used")

# Custom colors
observed_color <- "#437098"    # Dark blue
absent_color <- "white"        # Light beige
absent_outline <- "#D3D3D3"    # Light gray for outline
used_color <- "#E89275"        # Red for used years

# Plot the data
ggplot(penguin_abundance_complete, aes(x = year, y = site_id)) +
  geom_tile(aes(fill = presence_indicator, color = presence_indicator), linewidth = 0.2) +
  scale_fill_manual(values = c("Observed" = observed_color, "Absent" = absent_color), 
                    name = "Presence") +
  scale_color_manual(values = c("Observed" = observed_color, "Absent" = absent_outline), 
                     name = "Presence", guide = "none") +
  geom_tile(data = penguin_abundance_used, 
            aes(x = year, y = site_id, color = "Used in Model"), linewidth = 0.5, fill = NA) +
  scale_color_manual(values = c("Used in Model" = used_color), name = "Used in Model") +
  theme_minimal(base_size = 14) +
  labs(title = "Consecutive Sites Observed and Used in the Study",
       x = "Year",
       y = "Site ID") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10, margin = margin(t = 5, r = 5)), # Added margin for spacing
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )
