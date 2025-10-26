# ============================================================================
# data_processing_helpers.R
# Data utilities for penguin-sea ice analysis
#
# Required packages: dplyr, lubridate
# ============================================================================

#' Identify latitude column in dataframe
#'
#' @param df Data frame to search
#' @return String name of latitude column, or NULL if not found
pick_lat_col <- function(df) {
  cand <- c("lat", "latitude", "colony_lat", "site_lat", "lat_dd")
  hit <- cand[cand %in% names(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

#' Safe year extraction from date strings
#'
#' @param x Character vector of dates
#' @return Numeric vector of years
safe_year_from_date <- function(x) {
  x <- as.character(x)
  y <- suppressWarnings(year(mdy(x)))
  idx <- is.na(y)
  if (any(idx)) y[idx] <- suppressWarnings(year(dmy(x[idx])))
  idx <- is.na(y)
  if (any(idx)) y[idx] <- suppressWarnings(year(ymd(x[idx])))
  y
}

#' Aggregate metrics to yearly values
#'
#' Handles both daily and monthly input data
#'
#' @param df Data frame with metric data
#' @param metric_col String name of metric column to aggregate
#' @return Data frame with year and metric_value columns
yearly_aggregate_metric <- function(df, metric_col) {
  if ("date" %in% names(df)) {
    # Daily data → aggregate to year
    out <- df %>%
      mutate(year = safe_year_from_date(date)) %>%
      group_by(year) %>%
      summarise(metric_value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  } else if ("year" %in% names(df)) {
    # Already aggregated (monthly/seasonal) → average to year
    out <- df %>%
      group_by(year) %>%
      summarise(metric_value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  } else {
    stop("Input lacks both 'date' and 'year' columns; cannot aggregate.")
  }
  out
}