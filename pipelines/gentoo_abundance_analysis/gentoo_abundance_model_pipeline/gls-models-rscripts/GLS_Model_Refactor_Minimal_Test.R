# ================================================================
# Regional GLS tests across metrics, home-range sizes, thresholds,
# and lags (1..5), including persistence and SD metrics.
# ================================================================

library(tidyverse)
library(data.table)
library(nlme)
library(lubridate)

# -------------------------
# Helpers
# -------------------------
.pick_lat_col <- function(df) {
  cand <- c("lat", "latitude", "colony_lat", "site_lat", "lat_dd")
  hit  <- cand[cand %in% names(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

.safe_year_from_date <- function(x) {
  x <- as.character(x)
  y <- suppressWarnings(year(mdy(x)))
  idx <- is.na(y)
  if (any(idx)) y[idx] <- suppressWarnings(year(dmy(x[idx])))
  idx <- is.na(y)
  if (any(idx)) y[idx] <- suppressWarnings(year(ymd(x[idx])))
  y
}

# Given a data.frame with either daily "date" or columns "year" (and maybe "month"),
# return a yearly summary with column `metric_value`
.yearly_aggregate_metric <- function(df, metric_col) {
  if ("date" %in% names(df)) {
    # Daily → year
    out <- df %>%
      mutate(year = .safe_year_from_date(date)) %>%
      group_by(year) %>%
      summarise(metric_value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  } else if ("year" %in% names(df)) {
    # Monthly/seasonal rows → average to year
    # (If "month" exists with 6..9 & "Season-wide", this still averages across rows in that year.)
    out <- df %>%
      group_by(year) %>%
      summarise(metric_value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  } else {
    stop("Input lacks both 'date' and 'year' columns; cannot aggregate.")
  }
  out
}

# -------------------------
# Main function
# -------------------------
test_regional_gls <- function(penguin_data_path,
                              metric_files,
                              output_dir = NULL,
                              central_cutoff_deg = 63.2,
                              lags = 1:5,
                              min_rows_for_fit = 30,
                              metrics_priority = c(
                                # daily SIC/extent file
                                "mean_sic","sd_sic","ice_extent_km2",
                                # duration/persistence file
                                "mean_duration","sd_duration",
                                "mean_persistence","sd_persistence"
                              )) {
  
  # ---- Load penguins & classify regions
  penguins_raw <- fread(penguin_data_path)
  lat_col <- .pick_lat_col(penguins_raw)
  if (is.null(lat_col)) {
    stop("No latitude column found in penguin data. Expected one of: lat, latitude, colony_lat, site_lat, lat_dd.")
  }
  
  penguins <- penguins_raw %>%
    mutate(
      year = 1970 + season - 1,
      .lat_raw = .data[[lat_col]],
      .lat_sgn = if_else(.lat_raw < 0, .lat_raw, -abs(.lat_raw)),
      region   = if_else(.lat_sgn <= -central_cutoff_deg, "Central_WAP", "Bransfield")
    ) %>%
    filter(growth_rate <= 3, !is.na(growth_rate))
  
  # ---- Collect all results
  res_list <- list()
  
  for (f in metric_files) {
    dat <- fread(f)
    
    # Normalize presence of grouping columns
    if (!"HomeRangeSize" %in% names(dat)) dat[, HomeRangeSize := "Unknown"]
    if (!"Threshold"     %in% names(dat)) dat[, Threshold := NA_real_]
    
    # Which metric columns are in THIS file?
    metric_cols <- intersect(metrics_priority, names(dat))
    if (length(metric_cols) == 0) {
      warning(sprintf("No recognized metric columns in file: %s", f))
      next
    }
    
    home_sizes <- sort(unique(dat$HomeRangeSize))
    thresholds <- sort(unique(dat$Threshold))
    
    for (h in home_sizes) {
      for (th in thresholds) {
        
        df_ht <- dat %>% filter(HomeRangeSize == h, is.na(Threshold) | Threshold == th)
        if (nrow(df_ht) == 0) next
        
        for (metric_col in metric_cols) {
          
          # Map metric_col to a friendly name
          metric_name <- dplyr::case_when(
            metric_col == "mean_sic"         ~ "SIC",
            metric_col == "sd_sic"           ~ "SIC_SD",
            metric_col == "ice_extent_km2"   ~ "Extent",
            metric_col == "mean_duration"    ~ "Duration",
            metric_col == "sd_duration"      ~ "Duration_SD",
            metric_col == "mean_persistence" ~ "Persistence",
            metric_col == "sd_persistence"   ~ "Persistence_SD",
            TRUE                              ~ metric_col
          )
          
          # Build yearly time series for this metric
          d_year <- tryCatch(
            .yearly_aggregate_metric(df_ht, metric_col),
            error = function(e) {
              warning(sprintf("Skipping %s (HomeRange=%s, Thresh=%s, Metric=%s): %s",
                              basename(f), h, as.character(th), metric_col, e$message))
              return(NULL)
            }
          )
          if (is.null(d_year)) next
          
          # Evaluate across lags
          for (lag_y in lags) {
            
            d_lag <- d_year %>%
              mutate(year = year + lag_y)
            
            test_data <- penguins %>%
              left_join(d_lag, by = "year") %>%
              filter(!is.na(metric_value))
            
            if (nrow(test_data) < min_rows_for_fit) {
              # too few rows to fit reliably
              next
            }
            
            # ML fits for valid LRTs
            m1 <- gls(growth_rate ~ metric_value,
                      data = test_data,
                      correlation = corAR1(form = ~1 | site_id),
                      method = "ML")
            
            m2 <- gls(growth_rate ~ metric_value + region,
                      data = test_data,
                      correlation = corAR1(form = ~1 | site_id),
                      method = "ML")
            
            m3 <- gls(growth_rate ~ metric_value * region,
                      data = test_data,
                      correlation = corAR1(form = ~1 | site_id),
                      method = "ML")
            
            lr_main <- anova(m1, m2)
            lr_int  <- anova(m2, m3)
            
            tt <- summary(m1)$tTable
            pooled_slope   <- unname(tt["metric_value", "Value"])
            pooled_slope_p <- unname(tt["metric_value", "p-value"])
            
            yrs <- range(test_data$year, na.rm = TRUE)
            
            res_list[[length(res_list) + 1]] <- tibble(
              File                = basename(f),
              Metric              = metric_name,
              MetricColumn        = metric_col,
              HomeRange           = h,
              Threshold           = th,
              Lag                 = lag_y,
              N_Obs               = nrow(test_data),
              N_Colonies          = n_distinct(test_data$site_id),
              Year_Start          = yrs[1],
              Year_End            = yrs[2],
              MainEffect_p        = lr_main$`p-value`[2],
              Interaction_p       = lr_int$`p-value`[2],
              Pooled_Slope        = pooled_slope,
              Pooled_Slope_p      = pooled_slope_p,
              AIC_Combined        = AIC(m1),
              AIC_Region          = AIC(m2),
              AIC_Interaction     = AIC(m3),
              LogLik_Combined     = as.numeric(logLik(m1)),
              LogLik_Region       = as.numeric(logLik(m2)),
              LogLik_Interaction  = as.numeric(logLik(m3))
            )
          } # end lag
        }   # end metric_col
      }     # end threshold
    }       # end home size
  }         # end files loop
  
  out_df <- dplyr::bind_rows(res_list)
  
  if (!is.null(output_dir) && nrow(out_df) > 0) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    out_path <- file.path(output_dir, "Regional_Test_Summary_AllMetrics_AllLags.csv")
    readr::write_csv(out_df, out_path)
    message("Saved: ", out_path)
  }
  
  out_df
}

# -----------------------
# Example usage
# -----------------------
penguin_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ_temporary/gentoo-abundance-model/modeled_gentoo_parameters.csv"

metric_files <- c(
  "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ_temporary/gentoo-abundance-model/metric-calculation-csv/daily_sic_statistics.csv",
  "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ_temporary/gentoo-abundance-model/metric-calculation-csv/sea_ice_duration_persistence_stats.csv"
)

output_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ_temporary/RStudioProject/Results"

res_all <- test_regional_gls(
  penguin_data_path = penguin_path,
  metric_files      = metric_files,
  output_dir        = output_path,
  central_cutoff_deg = 63.2,
  lags               = 1:5,      # <-- evaluate 1..5 years prior
  min_rows_for_fit   = 30
)
