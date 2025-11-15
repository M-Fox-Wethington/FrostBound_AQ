
###############################################################################
# Load Necessary Libraries
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
library(mgcv)       # for GAMM models

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

# load(file.path(data_dir, "soi_data_long.rds"))

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
# Example Diagnostic Plots for an nlme Model (Random Slopes)
###############################################################################
# Suppose we have a random-slope model named model_lme_slopes
par(mfrow = c(2, 2))
plot(model_lme_slopes)         # Residuals vs. Fitted, etc.
qqnorm(resid(model_lme_slopes))
qqline(resid(model_lme_slopes))

# Residual vs. Fitted and QQ-Plot (for your final random-intercept model)
par(mfrow=c(1,2))
plot(model_random_intercept, which = 1)  # Residuals vs. Fitted
qqnorm(resid(model_random_intercept, type = 'p'))  
qqline(resid(model_random_intercept, type = 'p'))  # QQ-Plot


###############################################################################
# Fit Mixed-Effects Models for Different Lag Terms using nlme
###############################################################################
model_list <- list()
for (lag in 1:5) {
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
  cat("Lag", lag, "AIC:", AIC(model), "\n")
}

###############################################################################
# Refit a Model Using lmer (for Bootstrapping with bootMer)
###############################################################################
model_lmer <- lmer(growth_rate ~ winter_soi_lag1 + (winter_soi_lag1 | site),
                   data = merged_data, REML = TRUE)
summary(model_lmer)  # May give "boundary (singular) fit" warnings

boot_results <- bootMer(
  model_lmer,
  FUN = function(mod) fixef(mod),
  nsim = 1000
)
print(boot_results)

###############################################################################
# Geographical Visualization of Winter SOI (lag1)
###############################################################################
ggplot(merged_data_complete, aes(x = longitude, y = latitude)) +
  geom_point(aes(color = winter_soi_lag1), size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Geographical Distribution of Winter SOI (lag1)",
       x = "Longitude", y = "Latitude",
       color = "Winter SOI (lag1)") +
  theme_minimal()

###############################################################################
# Fit a GAMM to Explore Nonlinear Effects
###############################################################################
gamm_model <- gamm(
  growth_rate ~ s(winter_soi_lag1),
  random = list(site = ~1),
  correlation = corAR1(form = ~ year | site),
  data = merged_data,
  na.action = na.omit
)
summary(gamm_model$gam)
plot(gamm_model$gam, pages = 1)

###############################################################################
# Compare Random-Intercept vs. Random-Slope Models
###############################################################################
# 1) Random-intercept-only model
model_random_intercept <- lme(
  fixed = growth_rate ~ winter_soi_lag1,
  random = ~ 1 | site,
  correlation = corAR1(form = ~ year | site),
  data = merged_data,
  na.action = na.omit,
  control = lmeControl(opt = "optim", msMaxIter = 500)
)
summary(model_random_intercept)
cat("Random-Intercept Model AIC:", AIC(model_random_intercept), "\n")

# 2) Compare via likelihood ratio test
anova_results <- anova(model_random_intercept, model_lme_slopes)
print(anova_results)


###############################################################################
# Visualize
###############################################################################


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

library(ggplot2)
ggplot(df_plot, aes(x = fitted, y = observed)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Observed vs. Fitted Growth Rates",
    x = "Fitted",
    y = "Observed"
  ) +
  theme_minimal()


# Convert fitted and residuals to data frame, preserving row names
fitted_df <- data.frame(fitted = fitted_vals, resids = resids)
fitted_df$observed <- fitted_df$fitted + fitted_df$resids
fitted_df$row_id   <- row.names(fitted_df)

# Convert your original 'merged_data' to have row names that match
# (This requires that merged_data was used directly by the model without reordering.)
model_data <- merged_data
row.names(model_data) <- as.character(seq_len(nrow(model_data)))

# Join them by row names
plot_df <- merge(model_data, fitted_df, by.x = "row.names", by.y = "row_id")

# Now you can plot 'observed' vs. 'fitted', colored by site or year, etc.
ggplot(plot_df, aes(x = fitted, y = observed, color = as.factor(site))) +
  geom_point(alpha = 0.5) +
  theme_minimal()


###############################################################################
# Caterpillar Plots for Random Intercepts
###############################################################################


# Extract random intercepts
ranef_int <- ranef(model_random_intercept)
ranef_df <- data.frame(site = rownames(ranef_int),
                       intercept = ranef_int[[1]])

# Order by intercept estimate
ranef_df <- ranef_df %>%
  arrange(intercept) %>%
  mutate(order = row_number())

ggplot(ranef_df, aes(x = intercept, y = reorder(site, intercept))) +
  geom_point() +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  labs(title="Caterpillar Plot: Random Intercepts by Colony",
       x="Random Intercept Estimate", y="Colony") +
  theme_minimal()


final_fixef <- as.data.frame(summary(model_random_intercept)$tTable)
# Then you can convert it into a nicely formatted table
library(knitr)
kable(final_fixef, digits = 4,
      caption = "Fixed Effects for the Final Random-Intercept Model")


###############################################################################
# End of Script
###############################################################################
