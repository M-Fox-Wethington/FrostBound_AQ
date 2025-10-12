# ── Coefficient‐Trend Line Plots for Gentoo Penguin Growth‐Rate Models ──────────

# 1. Load libraries
library(tidyverse)   # ggplot2, dplyr, purrr, readr, etc.

# 2. Read & combine all CSVs
data_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/gentoo-model-results/subset-for-tables/R-Visualization"
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

df <- files %>%
  set_names() %>%
  map_dfr(read_csv, .id = "source_file")

# 3. Tidy data and rename for plotting
df2 <- df %>%
  rename(
    metric    = Metric,
    hr_size   = `Home Range Size`,
    lag       = Lag,
    estimate  = Coefficient,
    ci_lower  = `CI Lower`,
    ci_upper  = `CI Upper`,
    p_value   = `Bonferroni Adjusted p-Value`
  ) %>%
  # Convert home‐range to factor for color grouping
  mutate(hr_size = factor(hr_size, levels = sort(unique(hr_size))))

# 4. Build the trend‐line plot
trend_plot <- ggplot(df2, aes(x = lag, y = estimate, group = hr_size, color = hr_size)) +
  # line connecting estimates across lags
  geom_line(size = 1) +
  # points at each lag
  geom_point(size = 2) +
  # vertical error bars for 95% CI
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.2) +
  # one panel per metric
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  # zero‐effect reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  # axis and legend labels
  labs(
    title = "Coefficient Trends over Lags by Home‑Range Size",
    x     = "Lag (years)",
    y     = "Coefficient Estimate (95% CI)",
    color = "Home‑Range Size (km)"
  ) +
  theme_bw() +
  theme(
    strip.background   = element_rect(fill = "grey90"),
    strip.text         = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

# 5. Display the plot
print(trend_plot)

# 6. (Optional) Save to file
ggsave(
  filename = "coef_trend_line_plots.png",
  plot     = trend_plot,
  width    = 8,
  height   = 10,
  dpi      = 300
)
