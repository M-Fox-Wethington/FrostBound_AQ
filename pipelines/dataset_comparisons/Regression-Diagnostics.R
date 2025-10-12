
# Make sure rasters are compatible
if (!compareGeom(amsr_monthly, nsidc_monthly)) {
  stop("Rasters are not aligned")
}

# Convert rasters to vectors
amsr_values <- values(amsr_monthly)
nsidc_values <- values(nsidc_monthly)

# Create a combined data frame (each row corresponds to a pixel across all layers)
combined_df <- data.frame(AMSR = as.vector(amsr_values), NSIDC = as.vector(nsidc_values))

# Remove NA values to ensure clean regression analysis
clean_df <- na.omit(combined_df)

# Perform linear regression using the clean data
model <- lm(AMSR ~ NSIDC, data = clean_df)
summary(model)


plot(model)


# Plot NSIDC vs AMSR values
plot(clean_df$NSIDC, clean_df$AMSR, main="Scatter Plot of AMSR vs NSIDC", xlab="NSIDC", ylab="AMSR",
     pch=19, col="blue")  # Ensure this runs without error

# Add regression line
abline(model, col="red")



library(e1071)  # for skewness and kurtosis functions

# Assuming 'model_wls' is your model
residuals <- resid(model)

summary_stats <- summary(residuals)
skewness_value <- skewness(residuals)
kurtosis_value <- kurtosis(residuals)

# Create a summary list
residual_summary <- list(
  Mean = mean(residuals),
  Median = median(residuals),
  Standard_Deviation = sd(residuals),
  Skewness = skewness_value,
  Kurtosis = kurtosis_value
)

library(car)  # for outlier tests and diagnostics

# Calculate Cook's distance to identify influential cases
cooks_distances <- cooks.distance(model)

# High leverage points
hat_values <- hatvalues(model)

# Combine into a summary
diagnostic_summary <- data.frame(
  Cooks_Distance = cooks_distances,
  Leverage = hat_values
)

library(car)  # for outlier tests and diagnostics

# Calculate Cook's distance to identify influential cases
cooks_distances <- cooks.distance(model)

# High leverage points
hat_values <- hatvalues(model)

# Combine into a summary
diagnostic_summary <- data.frame(
  Cooks_Distance = cooks_distances,
  Leverage = hat_values
)


library(ggplot2)

# Basic Residual Plot
ggplot(data = diagnostic_summary, aes(x = seq_along(Cooks_Distance), y = Cooks_Distance)) +
  geom_line() +
  labs(title = "Cooks Distance Plot", x = "Observation", y = "Cooks Distance")

# You can include threshold lines or specific annotations based on the diagnostics


full_summary <- list(
  Residual_Summary = residual_summary,
  Diagnostic_Summary = diagnostic_summary
)

# If you want to present this in a readable format
# print(full_summary)

n <- nrow(clean_df)  # Total number of observations
threshold <- 4 / n

num_influential_cases <- sum(cooks_distances > threshold)

proportion_influential <- num_influential_cases / n

diagnostic_summary$Threshold_Cooks_Distance <- threshold
diagnostic_summary$Num_Influential_Cases <- num_influential_cases
diagnostic_summary$Proportion_Influential <- proportion_influential

print(paste("Number of influential cases (Cook's Distance > ", threshold, "): ", num_influential_cases))
print(paste("Proportion of influential cases: ", proportion_influential))


# Assuming 'model' is the  fitted weighted linear regression model
fitted_values <- fitted(model)
residuals <- resid(model)
standardized_residuals <- rstandard(model)

# Create a data frame for plotting
diagnostic_df <- data.frame(Fitted=fitted_values, Residuals=residuals, Standardized=standardized_residuals)

# Plot Fitted vs. Residuals
ggplot(diagnostic_df, aes(x=Fitted, y=Residuals)) +
  geom_point() +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  labs(title="Fitted vs. Residuals", x="Fitted Values", y="Residuals") +
  theme_minimal()

# Plot Fitted vs. Standardized Residuals
ggplot(diagnostic_df, aes(x=Fitted, y=Standardized)) +
  geom_point() +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  labs(title="Fitted vs. Standardized Residuals", x="Fitted Values", y="Standardized Residuals") +
  theme_minimal()



