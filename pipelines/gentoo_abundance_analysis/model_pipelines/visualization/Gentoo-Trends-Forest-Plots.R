# ── Forest Plot for Main Metrics (no NA estimates) ──────────────────────────
forest_main <- main_df %>%
  filter(!is.na(estimate)) %>%                  # drop the empty combos
  ggplot(aes(x = estimate, y = hr_size, color = lag)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(
    xmin = estimate - 1.96 * StdError,
    xmax = estimate + 1.96 * StdError
  ),
  height = 0.2, position = position_dodge(width = 0.7),
  na.rm = TRUE) +                             # also suppress NA warnings
  geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
  facet_grid(metric ~ threshold, scales = "free_y") +
  scale_color_viridis_d(name = "Lag (years)", end = 0.8) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.background   = element_rect(fill = "grey90", color = NA),
    strip.text         = element_text(face = "bold"),
    legend.position    = "bottom"
  ) +
  labs(
    # title = "Forest Plot: Main Sea Ice Metrics",
    x     = "Coefficient (95% CI)",
    y     = "Home‑Range Size (km)"
  )

print(forest_main)
ggsave("forest_main_metrics.png", forest_main, width = 10, height = 8, dpi = 300)


# ── Forest Plot for Variability Metrics (no NA estimates) ──────────────────
forest_var <- var_df %>%
  filter(!is.na(estimate)) %>%                  # drop the empty combos
  ggplot(aes(x = estimate, y = hr_size, color = lag)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(
    xmin = estimate - 1.96 * StdError,
    xmax = estimate + 1.96 * StdError
  ),
  height = 0.2, position = position_dodge(width = 0.7),
  na.rm = TRUE) +                             # also suppress NA warnings
  geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
  facet_grid(metric ~ threshold, scales = "free_y") +
  scale_color_viridis_d(name = "Lag (years)", end = 0.8) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.background   = element_rect(fill = "grey90", color = NA),
    strip.text         = element_text(face = "bold"),
    legend.position    = "bottom"
  ) +
  labs(
    # title = "Forest Plot: Variability Metrics",
    x     = "Coefficient (95% CI)",
    y     = "Home‑Range Size (km)"
  )

print(forest_var)
ggsave("forest_variability_metrics.png", forest_var, width = 10, height = 8, dpi = 300)




# ── 9. Per‑Metric Forest Plots ────────────────────────────────────────────────

# ── Derive base‐metric names dynamically ────────────────────────────────────────
base_metrics <- df2 %>%
  pull(metric) %>% 
  # remove any " Variability" suffix
  str_remove(" Variability$") %>% 
  # keep them unique
  unique()

# e.g. base_metrics might be:
# [1] "Sea Ice Concentration" "Sea Ice Extent"        
#     "Open Water Frequency" "Duration"              

# ── Then loop as before ───────────────────────────────────────────────────────
for (bm in base_metrics) {
  df_sub <- df2 %>%
    filter(metric %in% c(bm, paste0(bm, " Variability"))) %>%
    filter(!is.na(estimate))
  
  p <- ggplot(df_sub, aes(x = estimate, y = hr_size, color = lag)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(
      xmin = estimate - 1.96 * StdError,
      xmax = estimate + 1.96 * StdError
    ),
    height = 0.2, position = position_dodge(width = 0.7),
    na.rm = TRUE) +
    geom_point(position = position_dodge(width = 0.7), size = 2, na.rm = TRUE) +
    facet_grid(metric ~ threshold, scales = "free_y", switch = "y") +
    scale_color_viridis_d(name = "Lag (years)", end = 0.8) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.background   = element_rect(fill = "grey90", color = NA),
      strip.text         = element_text(face = "bold"),
      axis.title.y.left  = element_blank(),
      axis.text.y.left   = element_text(margin = margin(r = 10)),
      legend.position    = "bottom"
    ) +
    labs(
      title = paste0("Forest Plot: ", bm, " & Variability"),
      x     = "Coefficient (95% CI)",
      y     = "Home‑Range Size (km)"
    )
  
  print(p)
  
  safe_name <- gsub("[^[:alnum:]_]", "_", bm)
  ggsave(
    filename = paste0("forest_", safe_name, ".png"),
    plot     = p,
    width    = 10,
    height   = 6,
    dpi      = 300
  )
}
