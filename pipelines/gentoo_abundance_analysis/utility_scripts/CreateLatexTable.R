library(readr)
library(dplyr)
library(stringr)
library(glue)
library(purrr)

clean_and_numeric <- function(df) {
  num_cols <- c("Home Range Size","Threshold","Lag","Coefficient",
                "StdError","t-Value","Bonferroni Adjusted p-Value")
  # CI Lower/Upper added if missing
  for(col in intersect(c(num_cols, "CI Lower","CI Upper"), names(df))) {
    df[[col]] <- str_replace_all(df[[col]], "E", "e") |> as.numeric()
  }
  if(!"CI Lower" %in% names(df)) {
    df <- df %>%
      mutate(
        `CI Lower` = Coefficient - 1.96 * StdError,
        `CI Upper` = Coefficient + 1.96 * StdError
      )
  }
  df
}

build_rows <- function(df) {
  tbl <- df %>%
    select(Metric, `Home Range Size`, Threshold, Lag,
           Coefficient, StdError, `t-Value`,
           `Bonferroni Adjusted p-Value`, `CI Lower`, `CI Upper`) %>%
    rename(
      HomeRangeSize = `Home Range Size`,
      tValue        = `t-Value`,
      PBonf         = `Bonferroni Adjusted p-Value`,
      CIL           = `CI Lower`,
      CIU           = `CI Upper`
    )
  
  pmap_chr(tbl, function(Metric, HomeRangeSize, Threshold, Lag,
                         Coefficient, StdError, tValue, PBonf, CIL, CIU) {
    glue(
      "{Metric} & {HomeRangeSize} & ",
      "{sprintf('%.2f', Threshold)} & {Lag} & ",
      "{sprintf('%.3e', Coefficient)} & {sprintf('%.2e', StdError)} & ",
      "{sprintf('%.4f', tValue)} & {sprintf('%.4f', PBonf)} & ",
      "{sprintf('%.2e', CIL)} & {sprintf('%.2e', CIU)} \\\\"
    )
  })
}

build_latex_for_file <- function(csv_path, out_dir = ".") {
  # read + clean + sort
  df <- read_csv(csv_path, show_col_types = FALSE) %>%
    clean_and_numeric() %>%
    arrange(Metric, `Home Range Size`, Lag)
  
  # build rows
  rows     <- build_rows(df)
  body_tex <- paste(rows, collapse = "\n")
  
  # use filename for caption/label
  base   <- tools::file_path_sans_ext(basename(csv_path))
  caption <- glue("{base} results (sorted by Metric, home‑range size, lag) with 95\\% confidence intervals.")
  label   <- glue("tab:{str_replace_all(tolower(base), '[^a-z0-9]+','_')}")
  
  tex <- glue(
    "\\begin{{table}}[htbp]
\\centering
\\begin{{adjustbox}}{{max width=\\textwidth}}
\\begin{{tabular}}{{>{{\\raggedright\\arraybackslash}}m{{4cm}} c c c c c c c c c}}
\\toprule
Metric & Home Range Size & Threshold & Lag & Coefficient & StdError & t-Value & p-Value & CI Lower & CI Upper \\\\
\\midrule
{body_tex}
\\bottomrule
\\end{{tabular}}
\\end{{adjustbox}}
\\caption{{{caption}}}
\\label{{{label}}}
\\end{{table}}"
  )
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, paste0(base, ".tex"))
  write_file(tex, out_file)
  message("✓ Wrote ", out_file)
}

build_all_tables <- function(folder, out_sub = "latex_tables") {
  results_dir <- folder
  csvs        <- list.files(results_dir, pattern="\\.csv$", full.names=TRUE)
  out_dir     <- file.path(results_dir, out_sub)
  walk(csvs, build_latex_for_file, out_dir = out_dir)
}

# ─── Usage ─────────────────────────────────────────────────────────────────
build_all_tables(
  folder = "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/gentoo-model-results/subset-for-tables"
)

# csv_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/gentoo-model-results/subset-for-tables/Extent.csv"
# 
# df <- read_csv(csv_path, show_col_types=FALSE)
# print(names(df))

