# ============================================================================
# model_helpers.R
# GLS model fitting functions for penguin-sea ice analysis
#
# Required packages: nlme, dplyr
# ============================================================================

#' Fit GLS models for individual lag years (1-5)
#'
#' Tests lagged effects of sea ice metrics on penguin growth rates
#' using GLS with AR1 correlation structure
#'
#' @param penguin_data Data frame with penguin growth rates and site_id
#' @param metric_yearly Data frame with year and metric_value columns
#' @param metric_name String name of metric for labeling
#' @return List of significant model results (Bonferroni p < 0.05)
fit_gls_models_indiv <- function(penguin_data, metric_yearly, metric_name) {
  results <- list()
  
  for (lag in 1:5) {
    # Lag the metric data and merge with penguin data
    metric_lagged <- metric_yearly %>%
      mutate(year = year + lag)
    
    penguin_lagged <- penguin_data %>%
      left_join(metric_lagged, by = "year") %>%
      filter(!is.na(metric_value), !is.na(growth_rate), growth_rate <= 3)
    
    if (nrow(penguin_lagged) < 3) next
    
    cat(sprintf("    Fitting GLS model for lag %d year(s)\n", lag))
    
    formula <- as.formula("growth_rate ~ metric_value")
    model <- gls(formula, correlation = corAR1(form = ~1 | site_id), data = penguin_lagged)
    summary_model <- summary(model)
    bonferroni_p_value <- p.adjust(summary_model$tTable[2, 4], method = "bonferroni", n = 5)
    
    if (bonferroni_p_value < 0.05) {
      results[[paste("Lag_Indiv", lag)]] <- list(
        AIC = AIC(model),
        BIC = BIC(model),
        coefficients = summary_model$tTable,
        p_value = summary_model$tTable[2, 4],
        bonferroni_p_value = bonferroni_p_value,
        model = model,
        penguin_data = penguin_lagged,
        lag = lag,
        metric_name = metric_name
      )
    }
  }
  
  return(results)
}

#' Fit GLS models on pooled data (all regions combined)
#'
#' Same as fit_gls_models_indiv but optimized for pooled analysis
#'
#' @param penguin_data Data frame with penguin growth rates and site_id
#' @param metric_yearly Data frame with year and metric_value columns
#' @param metric_name String name of metric for labeling
#' @param report_all Logical, if TRUE return all models, not just significant
#' @return List of model results (significant only by default)
fit_gls_models_pooled <- function(penguin_data, metric_yearly, metric_name, 
                                  report_all = FALSE) {
  results <- list()
  
  for (lag in 1:5) {
    # Lag the metric data and merge with penguin data
    metric_lagged <- metric_yearly %>%
      mutate(year = year + lag)
    
    penguin_lagged <- penguin_data %>%
      left_join(metric_lagged, by = "year") %>%
      filter(!is.na(metric_value), !is.na(growth_rate), growth_rate <= 3)
    
    if (nrow(penguin_lagged) < 3) next
    
    cat(sprintf("    Fitting pooled GLS model for lag %d year(s)\n", lag))
    
    formula <- as.formula("growth_rate ~ metric_value")
    model <- gls(formula, correlation = corAR1(form = ~1 | site_id), data = penguin_lagged)
    summary_model <- summary(model)
    bonferroni_p_value <- p.adjust(summary_model$tTable[2, 4], method = "bonferroni", n = 5)
    
    # Store if significant OR if report_all = TRUE
    if (bonferroni_p_value < 0.05 || report_all) {
      results[[paste("Lag_Pooled", lag)]] <- list(
        AIC = AIC(model),
        BIC = BIC(model),
        coefficients = summary_model$tTable,
        p_value = summary_model$tTable[2, 4],
        bonferroni_p_value = bonferroni_p_value,
        model = model,
        penguin_data = penguin_lagged,
        lag = lag,
        metric_name = metric_name,
        is_significant = bonferroni_p_value < 0.05
      )
    }
  }
  
  return(results)
}