# ============================================================================
# Generate Figure S2: Growth Rate Spaghetti Plot
# "Growth Rate Trends Across Sites (1970-2023)"
#
# Recreates the growth_rate_spaghetti_all_labeled.png figure.
# Input: modeled_gentoo_parameters.csv (JAGS model output)
# Output: Figure_S2_growth_rate_spaghetti.png / .pdf
# ============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

# -- CONFIGURATION -----------------------------------------------------------

base_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica"

# Input: modeled growth rates from JAGS state-space model
input_csv <- file.path(base_dir,
  "FrostBound_AQ/Datasets/gentoo-abundance-model/inputs/modeled_gentoo_parameters.csv")

# Output directory
output_dir <- file.path(base_dir,
  "FrostBound_AQ/Publication_Figures/supplementary")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -- LOAD DATA ---------------------------------------------------------------

cat("Loading data from:", input_csv, "\n")
dat <- read_csv(input_csv, show_col_types = FALSE)

cat(sprintf("Loaded %d rows, %d unique sites\n", nrow(dat), n_distinct(dat$site_id)))

# -- PREPARE DATA ------------------------------------------------------------

# Filter out NA growth rates and ensure year column exists
plot_data <- dat %>%
  filter(!is.na(growth_rate)) %>%
  mutate(year = if ("year" %in% names(.)) year else 1970 + season - 1)

n_sites <- n_distinct(plot_data$site_id)
cat(sprintf("Plotting %d sites with growth rate data\n", n_sites))

# -- IDENTIFY OUTLIER PEAKS FOR LABELING -------------------------------------

# Find the single highest growth rate per site, then take the top outliers
site_peaks <- plot_data %>%
  group_by(site_id) %>%
  slice_max(growth_rate, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(growth_rate))

# Label the top 4 outlier peaks (matching original: GALE, BROW, BISC, CUVE)
n_labels <- 4
label_data <- site_peaks %>%
  slice_head(n = n_labels)

cat("Outlier peaks to label:\n")
print(label_data %>% select(site_id, year, growth_rate))

# -- CREATE PLOT -------------------------------------------------------------

p <- ggplot(plot_data, aes(x = year, y = growth_rate, group = site_id)) +
  # Individual site trajectories (light, semi-transparent)
  geom_line(alpha = 0.25, color = "steelblue", linewidth = 0.4) +
  # Outlier peak labels
  geom_text_repel(
    data = label_data,
    aes(x = year, y = growth_rate, label = site_id),
    inherit.aes = FALSE,
    size = 3.5,
    color = "grey30",
    nudge_y = 0.5,
    segment.color = "grey60",
    segment.size = 0.3,
    box.padding = 0.4
  ) +
  # Formatting
  labs(
    title = paste0("Growth Rate Trends Across ", n_sites, " Sites (1970-2023)"),
    x = "Year",
    y = "Growth Rate (Current / Previous Season)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  )

# -- SAVE --------------------------------------------------------------------

# PNG (300 DPI)
ggsave(
  file.path(output_dir, "Figure_S2_growth_rate_spaghetti.png"),
  plot = p, width = 10, height = 6, dpi = 300
)

# PDF
ggsave(
  file.path(output_dir, "Figure_S2_growth_rate_spaghetti.pdf"),
  plot = p, width = 10, height = 6, device = cairo_pdf
)

cat("\nFigure S2 saved to:", output_dir, "\n")
cat("  - Figure_S2_growth_rate_spaghetti.png\n")
cat("  - Figure_S2_growth_rate_spaghetti.pdf\n")
