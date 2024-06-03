# Install and load necessary libraries
# if (!require("rjags")) install.packages("rjags", dependencies = TRUE)
# if (!require("coda")) install.packages("coda", dependencies = TRUE)
library(rjags)
library(coda)
library(dplyr)
library(tidyr)
library(readr)
library(parallel)

model_file_path <- "D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/jags_model.jags"

# Write the JAGS model script to a file
sink(model_file_path)
cat("
model {

# process variance prior
sigma ~ dunif(0, 1)
tau <- pow(sigma, -2)

# breeding success: chick conversion priors
for (i in 1:chicks) {
  alpha[i] ~ dbeta(a, b)
}

for (i in 1:n_sites) {
  lz[i, s[i]] ~ dnorm(0, .001)
}

# intercept prior
beta ~ dunif(-.5, .5)

# site effects
for (i in 1:n_sites) {
  eta[i] ~ dnorm(0, tau_site)
}
sigma_site ~ dunif(0, 1)
tau_site <- pow(sigma_site, -2)

# season effects
for (t in 1:n_seasons) {
  epsilon[t] ~ dnorm(0, tau_season)
}
sigma_season ~ dunif(0, 1)
tau_season <- pow(sigma_season, -2)

# observation models
# nests
for (i in 1:nests) {
  y_n[i] ~ dnorm(mu_y_n[i], precision_n[i])
  mu_y_n[i] <- lz[site_n[i], season_n[i]] - 1/(2 * precision_n[i])
  y_n_new[i] ~ dnorm(mu_y_n[i], precision_n[i])
  y_n_sq[i] <- pow((y_n[i] - mu_y_n[i]), 2)
  y_n_sq_new[i] <- pow((y_n_new[i] - mu_y_n[i]), 2)
}

# chicks
for (i in 1:chicks) {
  N[i] <- 2 * round(exp(lz[site_c[i], season_c[i]]))
  z_c[i] ~ dbin(alpha[i], N[i])
  lz_c[i] <- log(z_c[i])
  y_c[i] ~ dnorm(mu_y_c[i], precision_c[i])
  mu_y_c[i] <- lz_c[i] - 1/(2 * precision_c[i])
  y_c_new[i] ~ dnorm(mu_y_c[i], precision_c[i])
  y_c_sq[i] <- pow((y_c[i] - mu_y_c[i]), 2)
  y_c_sq_new[i] <- pow((y_c_new[i] - mu_y_c[i]), 2)
}

# process model
# growth rate    
for (i in 1:n_sites) {
  for (t in 1:n_seasons) {
    zr[i, t] <- beta + eta[i] + epsilon[t]
    lza[i, t] <- lz[i, t] * w[i, t]
  }
}  

# initial year abundance
for (i in 1:n_sites) {
  y_i[i] ~ dnorm(mu_y_i[i], precision_i[i])
  mu_y_i[i] <- lz[i, s[i]] - 1/(2 * precision_i[i])
  y_i_new[i] ~ dnorm(mu_y_i[i], precision_i[i])
  y_i_sq[i] <- pow((y_i[i] - mu_y_i[i]), 2)
  y_i_sq_new[i] <- pow((y_i_new[i] - mu_y_i[i]), 2)
  
}

# abundance: initial year + 1 through max season
for (i in 1:n_sites) {
  for (t in (s[i] + 1):n_seasons) {
    lz[i, t] ~ dnorm(mu_lz[i, t], tau)
    mu_lz[i, t] <- lz[i, t - 1] + zr[i, t] - 1/(2 * tau)
  }
}

# abundance: min season through intial year - 1
for (i in 1:n_sites) {
  for (t in 1:(s[i] - 1)) {
    lz[i, s[i] - t] ~ dnorm(mu_lz[i, s[i] - t], tau)
    mu_lz[i, s[i] - t] <- lz[i, s[i] - t + 1] - zr[i, s[i] - t + 1] - 1/(2 * tau)
  }
}
 
# posterior predictive checks
y_n_sqs <- sum(y_n_sq[])
y_n_sqs_new <- sum(y_n_sq_new[])
y_i_sqs <- sum(y_i_sq[])
y_i_sqs_new <- sum(y_i_sq_new[])
y_c_sqs <- sum(y_c_sq[])
y_c_sqs_new <- sum(y_c_sq_new[])

# derived quantities
for (i in 1:n_sites) {
  for (t in 2:n_seasons) {
    l_a[i, t - 1] <- exp(lz[i, t] - lz[i, t - 1])
    l_p[i, t - 1] <- exp(zr[i, t])
    lw_a[i, t - 1] <- ifelse(sum(w[i, (t-1):t]) == 2, l_a[i, t - 1], 1)
    lw_p[i, t - 1] <- ifelse(sum(w[i, (t-1):t]) == 2, l_p[i, t - 1], 1)
  }
}

for (i in 1:n_sites) {
  x[i, 1:n_seasons] <- ifelse(sum(w[i, 1:n_seasons]) > 1, w[i, 1:n_seasons], rep(1, n_seasons)) 
  gl_a[i] <- ifelse(sum(w[i, 1:n_seasons]) > 1, pow(prod(lw_a[i, ]), (1/(sum(x[i, 1:n_seasons]) - 1))), 0)
  gl_p[i] <- ifelse(sum(w[i, 1:n_seasons]) > 1, pow(prod(lw_p[i, ]), (1/(sum(x[i, 1:n_seasons]) - 1))), 0)
}

}", fill = TRUE)
sink()

# Prepare the data list for JAGS model
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
  b = b
)


random_inits <- function(model_data) {
  seed = runif(1, 1, 100000)
  beta <- runif(1, -.025, .025)
  sigma_site <- runif(1, .025, .05)
  sigma_season <- runif(1, .05, .1)
  sigma <- runif(1, .05, .1)
  chicks <- model_data$chicks
  n_sites <- model_data$n_sites
  n_seasons <- model_data$n_seasons
  s <- model_data$s
  y_c <- model_data$y_c
  y_i <- model_data$y_i
  site_c <- model_data$site_c
  season_c <- model_data$season_c
  a <- model_data$a
  b <- model_data$b
  eta <- rnorm(n_sites, 0, sigma_site)
  epsilon <- rnorm(n_seasons, 0, sigma_season)
  alpha <- rbeta(chicks, a, b)
  lz <- zr <- array(NA, dim = c(n_sites, n_seasons))
  for (i in 1:n_sites) {
    lz[i, s[i]] <- mean(y_i[i], na.rm = TRUE)
    for (t in (s[i] + 1):n_seasons) {
      zr[i, t] <- beta + eta[i] + epsilon[t]
      lz[i, t] <- rnorm(1, lz[i, (t - 1)] + zr[i, t] - sigma^2 / 2, sigma)
    }
    for (t in 1:(s[i] - 1)) {
      zr[i, (s[i] - t + 1)] <- beta + eta[i] + epsilon[(s[i] - t + 1)]
      lz[i, (s[i] - t)] <- rnorm(1, lz[i, (s[i] - t + 1)] - zr[i, (s[i] - t + 1)] - sigma^2 / 2, sigma)
    }
  }
  
  z_c <- N <- NA
  for (i in 1:chicks) {
    if (lz[site_c[i], season_c[i]] < 0) lz[site_c[i], season_c[i]] <- 0
    N[i] <- 2 * round(exp(lz[site_c[i], season_c[i]]))
    z_c[i] <- base::max(rbinom(1, prob = alpha[i], N[i]), 1)
  }
  
  return(list(
    sigma = sigma,
    sigma_site = sigma_site,
    sigma_season = sigma_season,
    beta = beta,
    eta = eta,
    epsilon = epsilon,
    alpha = alpha,
    lz = lz,
    z_c = z_c,
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = seed))
}

save(random_inits, file = "D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/random_inits.rda")
expect_error(random_inits(model_data), NA)


n.chains <- 6
n.adapt <- 3000
n.update <- 300000
n.iter <- 200000
thin <- 200
cl <- makeCluster(n.chains)

cvars <- c("model_data", "n.adapt", "n.update", "n.iter", "thin", "params", "random_inits")
params <- c("beta", "sigma", "sigma_site", "sigma_season", "alpha", "epsilon", "eta", "z_c", "lz", 
            "gl_a", "l_a", "y_i_new", "y_n_new", "y_c_new", "y_n_sqs", "y_n_sqs_new", "y_i_sqs_new",
            "y_i_sqs", "y_c_sqs", "y_c_sqs_new", "lz_c", "lza")

parallel::clusterExport(cl, cvars)

out <- clusterEvalQ(cl, {
  library(rjags)
  inits <- random_inits(model_data)
  jm = jags.model("D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/jags_model.jags", data = model_data, n.chains = 1, n.adapt = n.adapt, 
                  inits = inits)
  update(jm, n.iter = n.update)
  zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = thin)
  return(as.mcmc(zm))
})
stopCluster(cl)
model_data_rinits_output = mcmc.list(out)  
save(model_data_rinits_output, file = "D:/Manuscripts_localData/FrostBound_AQ/Results/gentoo-abundance-model/model_data_rinits_output.rda")


model_data_rinits_output <- read

MCMCsummary(model_data_rinits_output, params = c("beta", "sigma", "sigma_site", "sigma_season"), 
            HPD = TRUE, hpd_prob = .95, round = 3)

#Posterior Predictive Checks
params <- c("y_i_sqs", "y_i_sqs_new", "y_n_sqs", "y_n_sqs_new", "y_c_sqs", "y_c_sqs_new")
MCMCsummary(model_data_rinits_output, params = params, n.eff = FALSE, round = 3)