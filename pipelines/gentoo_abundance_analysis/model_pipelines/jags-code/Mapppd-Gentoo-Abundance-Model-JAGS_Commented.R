# Install and load necessary libraries
# Optional installation lines are commented out. Uncomment to install if needed.
# if (!require("rjags")) install.packages("rjags", dependencies = TRUE)
# if (!require("coda")) install.packages("coda", dependencies = TRUE)
library(rjags)   # JAGS interface for Bayesian modeling
library(coda)    # Tools for MCMC output analysis
library(dplyr)
library(tidyr)
library(readr)
library(parallel) # Enables parallel processing

# Define the file path for the JAGS model script
model_file_path <- "D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/jags_model.jags"

# Write the JAGS model code to a text file
sink(model_file_path)
cat("
model {

# Define prior for process variance (sigma) and precision (tau)
sigma ~ dunif(0, 1)         # Uniform prior for sigma, allowing variation between 0 and 1
tau <- pow(sigma, -2)       # Precision is inverse squared sigma

# Breeding success prior for chick abundance
for (i in 1:chicks) {
  alpha[i] ~ dbeta(a, b)    # Beta distribution for breeding success, based on input parameters a and b
}

# Initial abundance for each site
for (i in 1:n_sites) {
  lz[i, s[i]] ~ dnorm(0, .001)  # Normally distributed initial abundance with high variance
}

# Define a prior for the intercept (beta) to allow variability in population growth
beta ~ dunif(-.5, .5)

# Site-specific effects (eta) for each site
for (i in 1:n_sites) {
  eta[i] ~ dnorm(0, tau_site)
}
sigma_site ~ dunif(0, 1)       # Uniform prior for site-level variation
tau_site <- pow(sigma_site, -2)

# Seasonal effects (epsilon) for each season
for (t in 1:n_seasons) {
  epsilon[t] ~ dnorm(0, tau_season)
}
sigma_season ~ dunif(0, 1)     # Uniform prior for season-level variation
tau_season <- pow(sigma_season, -2)

# Observation model for nest data
for (i in 1:nests) {
  y_n[i] ~ dnorm(mu_y_n[i], precision_n[i])   # Observed nest counts with specified precision
  mu_y_n[i] <- lz[site_n[i], season_n[i]] - 1/(2 * precision_n[i])  # Mean adjusted for precision
  y_n_new[i] ~ dnorm(mu_y_n[i], precision_n[i])   # Simulated nest counts for predictive checks
  y_n_sq[i] <- pow((y_n[i] - mu_y_n[i]), 2)       # Squared residuals for observed nest counts
  y_n_sq_new[i] <- pow((y_n_new[i] - mu_y_n[i]), 2)  # Squared residuals for simulated nest counts
}

# Observation model for chick data
for (i in 1:chicks) {
  N[i] <- 2 * round(exp(lz[site_c[i], season_c[i]]))   # Expected population size for chicks
  z_c[i] ~ dbin(alpha[i], N[i])                        # Binomially distributed chick abundance
  lz_c[i] <- log(z_c[i])                               # Log-transformed abundance
  y_c[i] ~ dnorm(mu_y_c[i], precision_c[i])            # Observed chick counts with specified precision
  mu_y_c[i] <- lz_c[i] - 1/(2 * precision_c[i])        # Mean for chick data, adjusted for precision
  y_c_new[i] ~ dnorm(mu_y_c[i], precision_c[i])        # Simulated chick counts for predictive checks
  y_c_sq[i] <- pow((y_c[i] - mu_y_c[i]), 2)            # Squared residuals for observed chick counts
  y_c_sq_new[i] <- pow((y_c_new[i] - mu_y_c[i]), 2)    # Squared residuals for simulated chick counts
}

# Process model for population growth
for (i in 1:n_sites) {
  for (t in 1:n_seasons) {
    zr[i, t] <- beta + eta[i] + epsilon[t]   # Growth rate influenced by site and seasonal effects
    lza[i, t] <- lz[i, t] * w[i, t]          # Adjusted abundance accounting for site-season presence
  }
}  

# Initial year abundance data model
for (i in 1:n_sites) {
  y_i[i] ~ dnorm(mu_y_i[i], precision_i[i])   # Observed abundance with precision for initial season
  mu_y_i[i] <- lz[i, s[i]] - 1/(2 * precision_i[i])  # Mean adjusted for precision
  y_i_new[i] ~ dnorm(mu_y_i[i], precision_i[i])      # Simulated initial counts for predictive checks
  y_i_sq[i] <- pow((y_i[i] - mu_y_i[i]), 2)          # Squared residuals for observed initial counts
  y_i_sq_new[i] <- pow((y_i_new[i] - mu_y_i[i]), 2)  # Squared residuals for simulated initial counts
}

# Population abundance dynamics for subsequent seasons
for (i in 1:n_sites) {
  for (t in (s[i] + 1):n_seasons) {
    lz[i, t] ~ dnorm(mu_lz[i, t], tau)              # Abundance in each season based on prior year
    mu_lz[i, t] <- lz[i, t - 1] + zr[i, t] - 1/(2 * tau)  # Process model with adjusted precision
  }
}

# Population abundance dynamics for pre-initial seasons
for (i in 1:n_sites) {
  for (t in 1:(s[i] - 1)) {
    lz[i, s[i] - t] ~ dnorm(mu_lz[i, s[i] - t], tau)       # Prior abundance modeled on next season
    mu_lz[i, s[i] - t] <- lz[i, s[i] - t + 1] - zr[i, s[i] - t + 1] - 1/(2 * tau)
  }
}
 
# Posterior predictive checks
y_n_sqs <- sum(y_n_sq[])      # Sum of squared residuals for nests (observed)
y_n_sqs_new <- sum(y_n_sq_new[])  # Sum of squared residuals for nests (simulated)
y_i_sqs <- sum(y_i_sq[])      # Sum of squared residuals for initial season (observed)
y_i_sqs_new <- sum(y_i_sq_new[])  # Sum of squared residuals for initial season (simulated)
y_c_sqs <- sum(y_c_sq[])      # Sum of squared residuals for chicks (observed)
y_c_sqs_new <- sum(y_c_sq_new[])  # Sum of squared residuals for chicks (simulated)

# Derived quantities for growth rate calculation
for (i in 1:n_sites) {
  for (t in 2:n_seasons) {
    l_a[i, t - 1] <- exp(lz[i, t] - lz[i, t - 1])   # Annual growth based on log-abundance change
    l_p[i, t - 1] <- exp(zr[i, t])                  # Growth rate (zr) transformation
    lw_a[i, t - 1] <- ifelse(sum(w[i, (t-1):t]) == 2, l_a[i, t - 1], 1)   # Weighted growth rate
    lw_p[i, t - 1] <- ifelse(sum(w[i, (t-1):t]) == 2, l_p[i, t - 1], 1)   # Weighted persistence rate
  }
}

# Geometric mean calculations for site-specific growth and persistence
for (i in 1:n_sites) {
  x[i, 1:n_seasons] <- ifelse(sum(w[i, 1:n_seasons]) > 1, w[i, 1:n_seasons], rep(1, n_seasons)) 
  gl_a[i] <- ifelse(sum(w[i, 1:n_seasons]) > 1, pow(prod(lw_a[i, ]), (1/(sum(x[i, 1:n_seasons]) - 1))), 0)
  gl_p[i] <- ifelse(sum(w[i, 1:n_seasons]) > 1, pow(prod(lw_p[i, ]), (1/(sum(x[i, 1:n_seasons]) - 1))), 0)
}

}", fill = TRUE)
sink() # Complete writing of the model to file