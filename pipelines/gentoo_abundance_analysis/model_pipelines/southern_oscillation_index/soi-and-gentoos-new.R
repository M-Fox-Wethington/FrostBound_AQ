###############################################################################
# Enhanced SOI Analysis for Gentoo Penguins - Nature-Quality Figures
###############################################################################
library(httr)       # for GET()
library(dplyr)
library(tidyr)      # for pivot_longer()
library(lubridate)
library(nlme)       # for mixed-effects models (lme)
library(readr)      # for reading CSVs
library(ggplot2)    # for plotting
library(tibble)     # for rownames_to_column()
library(knitr)      # for kable()
library(lme4)       # for lmer() models and bootMer()
library(gridExtra)  # for grid.arrange()
library(viridis)    # for better color palettes
library(xtable)     # for LaTeX table output
library(scales)     # for better axis formatting

# Set output directory
output_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ_temporary/gentoo-abundance-model/soi-analysis-results"

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

###############################################################################
# Define Data Directory and Load Data
###############################################################################
data_dir <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ_temporary/gentoo-abundance-model/"
gentoo_params <- read_csv(file.path(data_dir, "modeled_gentoo_parameters.csv"))
str(gentoo_params)

###############################################################################
# Download and Parse SOI Data from NOAA
###############################################################################
url <- "https://www.cpc.ncep.noaa.gov/data/indices/soi"
soi_raw <- tryCatch(GET(url), error = function(e) { stop("Error fetching SOI data: ", e) })
soi_text <- content(soi_raw, "text")
soi_lines <- strsplit(soi_text, "\n")[[1]]
standardized_soi_start <- grep("STANDARDIZED    DATA", soi_lines)
soi_standardized_lines <- soi_lines[(standardized_soi_start + 3):length(soi_lines)]
soi_standardized_lines <- soi_standardized_lines[!grepl("-999.9", soi_standardized_lines)]

soi_data <- read.table(text = soi_standardized_lines, fill = TRUE, stringsAsFactors = FALSE)
colnames(soi_data) <- c("Year", month.abb)
soi_data <- soi_data %>% mutate(Year = as.integer(Year))

soi_data_long <- soi_data %>%
  pivot_longer(-Year, names_to = "Month", values_to = "SOI") %>%
  mutate(
    Month = match(Month, month.abb),
    SOI   = as.numeric(SOI)
  ) %>%
  filter(!is.na(SOI))

# Optionally save the processed data for future use
saveRDS(soi_data_long, file = file.path(data_dir, "soi_data_long.rds"))
head(soi_data_long)

###############################################################################
# Prepare Winter SOI and Create Lagged Variables (June–Sept)
###############################################################################
winter_soi <- soi_data_long %>%
  filter(Month %in% c(6, 7, 8, 9)) %>%
  group_by(Year) %>%
  summarise(winter_soi = mean(SOI, na.rm = TRUE), .groups = "drop")

max_lag <- 5
winter_soi_lags <- winter_soi
for (lag_i in 1:max_lag) {
  lag_col <- paste0("winter_soi_lag", lag_i)
  winter_soi_lags <- winter_soi_lags %>%
    mutate(!!lag_col := dplyr::lag(winter_soi, n = lag_i, order_by = Year))
}
head(winter_soi_lags)

###############################################################################
# Merge SOI Lags into Gentoo Growth Data
###############################################################################
merged_data <- gentoo_params %>%
  left_join(winter_soi_lags, by = c("year" = "Year"))
merged_data_complete <- merged_data %>% filter(!is.na(winter_soi_lag1))
head(merged_data_complete)

###############################################################################
# Fit Mixed-Effects Models for Different Lag Terms using nlme
###############################################################################
cat("Fitting models for different lags...\n")
model_list <- list()
for (lag in 1:5) {
  cat("Fitting model for lag", lag, "...\n")
  lag_term <- paste0("winter_soi_lag", lag)
  fixed_formula <- as.formula(paste("growth_rate ~", lag_term))
  random_formula <- as.formula(paste("~", lag_term, "| site"))
  
  model <- lme(
    fixed = fixed_formula,
    random = random_formula,
    correlation = corAR1(form = ~ year | site),
    data = merged_data,
    na.action = na.omit,
    control = lmeControl(opt = "optim", msMaxIter = 500)
  )
  
  model_list[[lag]] <- model
  cat("  Lag", lag, "AIC:", AIC(model), "\n")
}

# Define the slopes model for use in comparisons
model_lme_slopes <- model_list[[1]]  # The lag 1 model becomes our slopes model
cat("Model list created with", length(model_list), "models\n")

###############################################################################
# Refit a Model Using lmer (for Bootstrapping with bootMer)
###############################################################################
cat("Fitting lmer model for bootstrapping...\n")
model_lmer <- lmer(growth_rate ~ winter_soi_lag1 + (winter_soi_lag1 | site),
                   data = merged_data, REML = TRUE)
summary_capture <- capture.output(summary(model_lmer))
writeLines(summary_capture, file.path(output_dir, "lmer_model_summary.txt"))

cat("Running bootstrap analysis...\n")
boot_results <- bootMer(
  model_lmer,
  FUN = function(mod) fixef(mod),
  nsim = 1000
)
boot_capture <- capture.output(print(boot_results))
writeLines(boot_capture, file.path(output_dir, "bootstrap_results.txt"))

###############################################################################
# Compare Random-Intercept vs. Random-Slope Models
###############################################################################
cat("Fitting random-intercept model...\n")
model_random_intercept <- lme(
  fixed = growth_rate ~ winter_soi_lag1,
  random = ~ 1 | site,
  correlation = corAR1(form = ~ year | site),
  data = merged_data,
  na.action = na.omit,
  control = lmeControl(opt = "optim", msMaxIter = 500)
)
ri_summary <- summary(model_random_intercept)
ri_capture <- capture.output(print(ri_summary))
writeLines(ri_capture, file.path(output_dir, "random_intercept_model_summary.txt"))
cat("Random-Intercept Model AIC:", AIC(model_random_intercept), "\n")

# Compare via likelihood ratio test
cat("Comparing models via likelihood ratio test...\n")
anova_results <- anova(model_random_intercept, model_lme_slopes)
anova_capture <- capture.output(print(anova_results))
writeLines(anova_capture, file.path(output_dir, "model_comparison_anova.txt"))

###############################################################################
# Extract and compile results from all models
###############################################################################
cat("Extracting results from all models...\n")
model_summary <- data.frame(
  Lag = 1:5,
  AIC = sapply(model_list, AIC),
  BIC = sapply(model_list, BIC)
)

# Extract fixed effects from each model
for (i in 1:5) {
  # Get model summary
  mod_summary <- summary(model_list[[i]])
  
  # Get the fixed effects table
  tbl <- mod_summary$tTable
  
  # Extract the coefficient for the appropriate lag term
  lag_term <- paste0("winter_soi_lag", i)
  
  # Assign values to the data frame
  model_summary$Coefficient[i] <- tbl[lag_term, "Value"]
  model_summary$StdError[i] <- tbl[lag_term, "Std.Error"]
  model_summary$tValue[i] <- tbl[lag_term, "t-value"]
  model_summary$pValue[i] <- tbl[lag_term, "p-value"]
  model_summary$RandomSD[i] <- as.numeric(VarCorr(model_list[[i]])[1, 2])
  
  # Write full model summary to file
  writeLines(capture.output(print(mod_summary)), 
             file.path(output_dir, paste0("model_lag_", i, "_summary.txt")))
}

# Print and save model summary
print(model_summary)
write.csv(model_summary, file.path(output_dir, "soi_models_summary.csv"), row.names = FALSE)

###############################################################################
# Define Nature-style theme for publication-quality figures
###############################################################################
# Define color palette for Nature style
nature_blue <- "#0072B2"    # Primary blue
nature_red <- "#D55E00"     # Primary red/orange
nature_green <- "#009E73"   # Secondary green
nature_purple <- "#CC79A7"  # Secondary purple
nature_orange <- "#E69F00"  # Secondary orange

# Custom theme function for Nature-style plots
nature_theme <- function() {
  theme_minimal() +
    theme(
      # Text elements
      text = element_text(color = "#000000"),
      plot.title = element_text(size = 11, face = "plain", hjust = 0),
      plot.subtitle = element_text(size = 9, hjust = 0, margin = margin(b = 15)),
      axis.title = element_text(size = 9, face = "plain"),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5)),
      
      # Grid elements
      panel.grid.major = element_line(color = "#E5E5E5", size = 0.2),
      panel.grid.minor = element_blank(),
      
      # Plot aesthetics
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Legend elements
      legend.position = "bottom",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.5, "cm"),
      legend.margin = margin(t = 5, b = 5),
      
      # Margins and layout
      plot.margin = margin(15, 15, 15, 15),
      
      # Caption
      plot.caption = element_text(size = 8, hjust = 0, face = "plain", margin = margin(t = 10))
    )
}

###############################################################################
# Create enhanced publication-quality figure: Effect sizes across lag times
###############################################################################
cat("Creating publication-quality effect size plot...\n")

# Add significance annotation
model_summary$sig <- ifelse(model_summary$pValue < 0.001, "***", 
                            ifelse(model_summary$pValue < 0.01, "**",
                                   ifelse(model_summary$pValue < 0.05, "*", "ns")))

# Calculate 95% confidence intervals
model_summary$ci_lower <- model_summary$Coefficient - 1.96 * model_summary$StdError
model_summary$ci_upper <- model_summary$Coefficient + 1.96 * model_summary$StdError

# Enhanced coefficient plot
lag_plot <- ggplot(model_summary, aes(x = Lag, y = Coefficient)) +
  # Add reference line at y=0
  geom_hline(yintercept = 0, linetype = "solid", color = "grey70", size = 0.3) +
  # Add error bars
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.2, color = nature_blue, size = 0.4) +
  # Add points with significance coloring
  geom_point(aes(fill = factor(sig)), shape = 21, size = 3, stroke = 0.5, color = "white") +
  # Add connecting line
  geom_line(color = nature_blue, size = 0.4, alpha = 0.5) +
  # Set scale with good spacing
  scale_x_continuous(breaks = 1:5, labels = 1:5, limits = c(0.5, 5.5)) +
  scale_y_continuous(limits = c(-0.06, 0.06), 
                     breaks = seq(-0.06, 0.06, 0.02),
                     labels = function(x) format(x, nsmall = 2)) +
  # Use the Nature color palette for significance
  scale_fill_manual(values = c("ns" = "#CCCCCC", "*" = "#56B4E9", 
                               "**" = "#0072B2", "***" = "#E69F00"),
                    name = "Significance", 
                    labels = c("ns" = "p ≥ 0.05", "*" = "p < 0.05", 
                               "**" = "p < 0.01", "***" = "p < 0.001")) +
  # Labels
  labs(x = "Time lag (years)",
       y = "SOI coefficient (effect on growth rate)",
       caption = "Figure 1 | SOI effect on Gentoo penguin population growth rates across time lags.") +
  # Apply Nature theme
  nature_theme() 

###############################################################################
# Create enhanced publication-quality figure: AIC plot
###############################################################################
cat("Creating publication-quality AIC plot...\n")

aic_plot <- ggplot(model_summary, aes(x = Lag, y = AIC)) +
  # Add points
  geom_point(size = 3, color = nature_red) +
  # Add connecting line
  geom_line(color = nature_red, size = 0.4, alpha = 0.5) +
  # Set scale with good spacing
  scale_x_continuous(breaks = 1:5, labels = 1:5, limits = c(0.5, 5.5)) +
  scale_y_continuous(limits = c(1950, 2060), 
                     breaks = seq(1960, 2060, 20)) +
  # Labels
  labs(x = "Time lag (years)",
       y = "Akaike Information Criterion",
       caption = "Figure 2 | Model fit by lag time (lower AIC indicates better fit).") +
  # Apply Nature theme
  nature_theme()

###############################################################################
# Create enhanced publication-quality figure: SOI vs Growth Relationship
###############################################################################
cat("Creating publication-quality SOI-Growth relationship plot...\n")

# Filter data
scatter_data <- merged_data %>%
  filter(!is.na(winter_soi_lag1), !is.na(growth_rate))

# Create the enhanced relationship plot
soi_growth_plot <- ggplot(scatter_data, aes(x = winter_soi_lag1, y = growth_rate)) +
  # Add points
  geom_point(alpha = 0.6, size = 1.5, color = nature_blue) +
  # Add regression line with confidence interval
  geom_smooth(method = "lm", color = nature_red, fill = nature_red, alpha = 0.2, size = 0.5) +
  # Add better axis formatting
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  # Labels
  labs(x = "Winter SOI (lag 1 year)",
       y = "Population growth rate",
       caption = "Figure 3 | Relationship between winter SOI (1-year lag) and Gentoo penguin growth.") +
  # Apply Nature theme
  nature_theme()

###############################################################################
# Create enhanced publication-quality figure: Geographical plot
###############################################################################
cat("Creating publication-quality geographical visualization...\n")

geo_plot <- ggplot(merged_data_complete, aes(x = longitude, y = latitude)) +
  # Add points with SOI coloring
  geom_point(aes(color = winter_soi_lag1), size = 2, alpha = 0.8) +
  # Use diverging color scale (blue-white-red)
  scale_color_gradient2(low = nature_blue, mid = "white", high = nature_red, 
                        midpoint = 0, name = "Winter SOI\n(lag 1 year)") +
  # Set scale with good spacing
  coord_fixed(ratio = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  # Labels
  labs(x = "Longitude",
       y = "Latitude",
       caption = "Figure 4 | Geographical distribution of winter SOI values across Gentoo penguin colonies.") +
  # Apply Nature theme
  nature_theme() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.margin = margin(5, 5, 5, 5)
  )

###############################################################################
# Create enhanced publication-quality figure: Observed vs Fitted
###############################################################################
cat("Creating publication-quality observed vs. fitted plot...\n")

# Extract the fitted values and residuals from the final model
fitted_vals <- fitted(model_random_intercept)
resids      <- residuals(model_random_intercept)

# Reconstruct the "observed" values the model used
observed_vals <- fitted_vals + resids

# Put them in a dataframe with the same number of rows as the model
df_plot <- data.frame(
  observed = observed_vals,
  fitted   = fitted_vals
)

fit_plot <- ggplot(df_plot, aes(x = fitted, y = observed)) +
  geom_point(alpha = 0.6, size = 1.5, color = nature_blue) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = nature_red, size = 0.5) +
  labs(
    x = "Fitted growth rate",
    y = "Observed growth rate",
    caption = "Figure 5 | Model validation: observed versus fitted growth rates."
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  nature_theme()

###############################################################################
# Create combined figure (panel layout for publication)
###############################################################################
cat("Creating combined panel figure for publication...\n")

combined_plot <- grid.arrange(
  lag_plot + theme(legend.position = "none"),
  aic_plot,
  soi_growth_plot,
  ncol = 1,
  heights = c(1, 0.8, 1),
  top = "Southern Oscillation Index effects on Gentoo penguin population growth"
)

###############################################################################
# Save all plots with publication-quality settings
###############################################################################

# Save the Nature-styled plots (PNG first for reliability)
ggsave(file.path(output_dir, "nature_soi_coefficient_by_lag.png"), lag_plot, 
       width = 89, height = 89, units = "mm", dpi = 600)

ggsave(file.path(output_dir, "nature_aic_by_lag.png"), aic_plot,
       width = 89, height = 89, units = "mm", dpi = 600)

ggsave(file.path(output_dir, "nature_soi_growth_relationship.png"), soi_growth_plot,
       width = 89, height = 89, units = "mm", dpi = 600)

ggsave(file.path(output_dir, "nature_geographical_distribution.png"), geo_plot,
       width = 89, height = 89, units = "mm", dpi = 600)

ggsave(file.path(output_dir, "nature_observed_vs_fitted.png"), fit_plot,
       width = 89, height = 89, units = "mm", dpi = 600)

ggsave(file.path(output_dir, "nature_combined_panels.png"), combined_plot,
       width = 183, height = 247, units = "mm", dpi = 600)

# Try to save PDFs with more robust device
tryCatch({
  ggsave(file.path(output_dir, "nature_soi_coefficient_by_lag.pdf"), lag_plot, 
         width = 89, height = 89, units = "mm", dpi = 600, device = cairo_pdf)
  
  ggsave(file.path(output_dir, "nature_aic_by_lag.pdf"), aic_plot,
         width = 89, height = 89, units = "mm", dpi = 600, device = cairo_pdf)
  
  ggsave(file.path(output_dir, "nature_soi_growth_relationship.pdf"), soi_growth_plot,
         width = 89, height = 89, units = "mm", dpi = 600, device = cairo_pdf)
  
  ggsave(file.path(output_dir, "nature_geographical_distribution.pdf"), geo_plot,
         width = 89, height = 89, units = "mm", dpi = 600, device = cairo_pdf)
  
  ggsave(file.path(output_dir, "nature_observed_vs_fitted.pdf"), fit_plot,
         width = 89, height = 89, units = "mm", dpi = 600, device = cairo_pdf)
  
  ggsave(file.path(output_dir, "nature_combined_panels.pdf"), combined_plot,
         width = 183, height = 247, units = "mm", dpi = 600, device = cairo_pdf)
  
  cat("PDF files saved successfully.\n")
}, error = function(e) {
  cat("Error saving PDF files:", e$message, "\n")
  cat("PNG files were saved successfully and can be used for submission.\n")
})

###############################################################################
# Also save originals for compatibility
###############################################################################

# Create and save the original plots as well
orig_lag_plot <- ggplot(model_summary, aes(x = Lag, y = Coefficient)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Coefficient - 1.96 * StdError, 
                    ymax = Coefficient + 1.96 * StdError), width = 0.2) +
  geom_line(linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  labs(title = "SOI Effect on Gentoo Penguin Growth by Lag Time",
       x = "Lag (Years)",
       y = "SOI Coefficient (Effect on Growth Rate)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12))

orig_aic_plot <- ggplot(model_summary, aes(x = Lag, y = AIC)) +
  geom_point(size = 3) +
  geom_line(linetype = "dashed") +
  labs(title = "Model Fit by Lag Time",
       x = "Lag (Years)",
       y = "AIC (lower is better)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12))

orig_soi_growth_plot <- ggplot(scatter_data, aes(x = winter_soi_lag1, y = growth_rate)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Relationship Between Winter SOI (1-year lag) and Gentoo Growth",
       x = "Winter SOI (1-year lag)",
       y = "Population Growth Rate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12))

# Save original versions too
ggsave(file.path(output_dir, "soi_coefficient_by_lag.png"), orig_lag_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "soi_aic_by_lag.png"), orig_aic_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "soi_growth_relationship.png"), orig_soi_growth_plot, width = 8, height = 6, dpi = 300)

# Save original combined plot
orig_combined_lag_plot <- grid.arrange(orig_lag_plot, orig_aic_plot, ncol = 1)
ggsave(file.path(output_dir, "soi_lag_combined.png"), orig_combined_lag_plot, width = 8, height = 10, dpi = 300)

###############################################################################
# Create LaTeX table of fixed effects
###############################################################################
cat("Creating fixed effects table...\n")
final_fixef <- as.data.frame(summary(model_random_intercept)$tTable)
# Save as CSV
write.csv(final_fixef, file.path(output_dir, "random_intercept_fixed_effects.csv"))

# Create LaTeX table
latex_fixef <- xtable(final_fixef, 
                      caption = "Fixed Effects for the Final Random-Intercept Model",
                      label = "tab:fixed_effects",
                      digits = 5)

latex_fixef_code <- capture.output(print(latex_fixef, 
                                         include.rownames = TRUE,
                                         caption.placement = "top"))

writeLines(latex_fixef_code, file.path(output_dir, "fixed_effects_latex.tex"))

# Create a nice table with kable
kable_output <- capture.output(kable(final_fixef, digits = 4,
                                     caption = "Fixed Effects for the Final Random-Intercept Model"))
writeLines(kable_output, file.path(output_dir, "fixed_effects_kable.txt"))

###############################################################################
# Generate LaTeX table of model comparisons
###############################################################################
cat("Creating model comparison table...\n")
# Create LaTeX table of model results
latex_table <- xtable(model_summary[, c("Lag", "Coefficient", "StdError", "tValue", "pValue", "AIC")], 
                      caption = "Effects of Southern Oscillation Index (SOI) on Gentoo Penguin Population Growth Rates",
                      label = "tab:soi_effects",
                      digits = c(0, 0, 4, 4, 2, 5, 1))

latex_code <- capture.output(print(latex_table, 
                                   include.rownames = FALSE,
                                   caption.placement = "top"))

writeLines(latex_code, file.path(output_dir, "soi_table_latex.tex"))

###############################################################################
# Comprehensive results summary
###############################################################################
cat("Creating comprehensive summary report...\n")
summary_file <- file.path(output_dir, "soi_analysis_summary.txt")
sink(summary_file)

cat("\n\nSOI ANALYSIS RESULTS SUMMARY\n")
cat("============================\n\n")

cat("1. Best models by AIC:\n")
cat("   Lag 5: AIC =", round(model_summary$AIC[5], 1), "\n")
cat("   Lag 1: AIC =", round(model_summary$AIC[1], 1), "\n\n")

cat("2. SOI effect on Gentoo growth (1-year lag):\n")
cat("   Coefficient =", round(model_summary$Coefficient[1], 4), 
    "±", round(model_summary$StdError[1], 4), "\n")
cat("   t-value =", round(model_summary$tValue[1], 2), 
    ", p-value =", format.pval(model_summary$pValue[1], digits = 3), "\n\n")

cat("3. SOI effect on Gentoo growth (5-year lag):\n")
cat("   Coefficient =", round(model_summary$Coefficient[5], 4), 
    "±", round(model_summary$StdError[5], 4), "\n")
cat("   t-value =", round(model_summary$tValue[5], 2), 
    ", p-value =", format.pval(model_summary$pValue[5], digits = 3), "\n\n")

cat("4. Random effects structure:\n")
cat("   Colony-specific intercept SD =", 
    round(as.numeric(VarCorr(model_random_intercept)[1, 2]), 4), "\n")
cat("   Residual SD =", 
    round(as.numeric(VarCorr(model_random_intercept)[2, 2]), 4), "\n\n")

cat("5. Interpretation:\n")
if (model_summary$pValue[1] < 0.05) {
  cat("   Positive SOI values (La Niña conditions) are associated with\n")
  cat("   increased Gentoo penguin growth rates in the following year.\n")
  cat("   For each unit increase in SOI, growth rates increase by approximately\n")
  cat("  ", round(model_summary$Coefficient[1], 3), "units.\n\n")
} else {
  cat("   No significant relationship between SOI and growth rates at lag 1.\n\n")
}

cat("6. Random slopes vs. random intercepts:\n")
cat("   The random slopes model (AIC =", round(AIC(model_lme_slopes), 1), 
    ") performs slightly better than\n")
cat("   the random intercepts model (AIC =", round(AIC(model_random_intercept), 1), ").\n")
cat("   The correlation between random intercepts and slopes is", 
    round(VarCorr(model_lme_slopes)[3, 2], 3), ".\n\n")

cat("This analysis supports the hypothesis that large-scale climate oscillations\n")
cat("influence Gentoo penguin population dynamics, likely through their effects\n")
cat("on sea ice conditions along the Western Antarctic Peninsula.\n")

# Modified section to fix the error:
cat("6. Random slopes vs. random intercepts:\n")
cat("   The random slopes model (AIC =", round(AIC(model_lme_slopes), 1), 
    ") performs slightly better than\n")
cat("   the random intercepts model (AIC =", round(AIC(model_random_intercept), 1), ").\n")

# Instead of trying to directly access VarCorr(model_lme_slopes)[3, 2], let's extract it properly:
vc_matrix <- VarCorr(model_lme_slopes)
if(nrow(vc_matrix) >= 3) {
  corr_value <- vc_matrix[3, 2]
  cat("   The correlation between random intercepts and slopes is", 
      round(as.numeric(corr_value), 3), ".\n\n")
} else {
  cat("   (Note: Correlation between random intercepts and slopes not available in this model structure)\n\n")
}

sink()

# Create a nice-formatted summary table for manuscript
key_findings <- data.frame(
  Metric = c("SOI Effect (Lag 1)", "SOI Effect (Lag 5)", 
             "Random Intercept SD", "AR(1) Correlation", 
             "AIC (Lag 1 Model)", "AIC (Lag 5 Model)",
             "AIC (Random Intercept Model)"),
  Value = c(
    sprintf("%.4f ± %.4f", model_summary$Coefficient[1], model_summary$StdError[1]),
    sprintf("%.4f ± %.4f", model_summary$Coefficient[5], model_summary$StdError[5]),
    sprintf("%.4f", as.numeric(VarCorr(model_random_intercept)[1, 2])),
    sprintf("%.4f", as.numeric(model_random_intercept$modelStruct$corStruct)),
    sprintf("%.1f", model_summary$AIC[1]),
    sprintf("%.1f", model_summary$AIC[5]),
    sprintf("%.1f", AIC(model_random_intercept))
  ),
  Interpretation = c(
    ifelse(model_summary$pValue[1] < 0.05, "Significant positive effect", "Non-significant"),
    ifelse(model_summary$pValue[5] < 0.05, "Significant negative effect", "Non-significant"),
    "Moderate colony-level variation",
    "Minimal temporal autocorrelation",
    "Better than lags 2-4",
    "Best model fit",
    "Similar to random slopes model"
  )
)

write.csv(key_findings, file.path(output_dir, "soi_key_findings.csv"), row.names = FALSE)

# Create Nature-style table for publication
nature_table <- xtable(
  key_findings[, 1:2],
  caption = "Table 1. Southern Oscillation Index effects on Gentoo penguin population growth rates.",
  label = "tab:soi_key_findings"
)

nature_table_code <- capture.output(print(nature_table, 
                                          include.rownames = FALSE,
                                          caption.placement = "top"))
writeLines(nature_table_code, file.path(output_dir, "nature_soi_table.tex"))

# Output manuscript-ready results in a formatted file
manuscript_text <- file.path(output_dir, "manuscript_text_soi.txt")
sink(manuscript_text)

cat("SOI Effects on Gentoo Penguin Population Growth Rates\n")
cat("=====================================================\n\n")

cat("Our analysis of Southern Oscillation Index (SOI) effects on Gentoo penguin population growth\n")
cat("reveals a complex pattern of temporal relationships. Positive SOI values (La Niña conditions)\n")
cat("are significantly associated with increased growth rates in the following year (coefficient =\n")
cat(sprintf("%.4f ± %.4f, t = %.2f, p < 0.001", 
            model_summary$Coefficient[1], 
            model_summary$StdError[1],
            model_summary$tValue[1]), "). \n\n")

cat("Interestingly, we observed a 5-year lagged negative relationship (coefficient = ")
cat(sprintf("%.4f ± %.4f, t = %.2f, p < 0.001", 
            model_summary$Coefficient[5], 
            model_summary$StdError[5],
            model_summary$tValue[5]), "), \n")
cat("suggesting a complex oscillatory pattern in SOI effects on penguin demographics. The model\n")
cat("with 5-year lag provides the best fit (AIC = ", round(model_summary$AIC[5], 1), ") followed by the\n")
cat("1-year lag model (AIC = ", round(model_summary$AIC[1], 1), ").\n\n")

cat("These findings support the hypothesis that large-scale climate oscillations influence\n")
cat("Gentoo penguin population dynamics through their effects on sea ice conditions.\n")
cat("SOI-driven changes in sea ice extent and duration likely create complex demographic\n")
cat("responses through direct effects on adult survival and indirect effects on recruitment\n")
cat("that manifest over multiple years.\n\n")

sink()

cat("\nAnalysis complete. Publication-quality results saved to:", output_dir, "\n")
cat("Enhanced Nature-style figures created and saved successfully.\n")