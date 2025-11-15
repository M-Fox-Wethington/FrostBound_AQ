# ============================================================================
# Regional Pooling Justification Analysis
# Quantifies regional differences to support pooling decision
# ============================================================================

library(dplyr)
library(data.table)

#' Analyze regional differences and generate pooling justification
#'
#' @param regional_comparison_path Path to Regional_Comparison_Summary.csv
#' @param interaction_results_path Path to interaction_tests_all.csv
#' @param output_dir Directory to save summary reports
#' @return List with summary statistics and decision justification
analyze_pooling_justification <- function(regional_comparison_path,
                                          interaction_results_path = NULL,
                                          output_dir = NULL) {
  
  cat("\n", rep("=", 70), "\n")
  cat("REGIONAL POOLING JUSTIFICATION ANALYSIS\n")
  cat(rep("=", 70), "\n\n")
  
  # Load regional comparison data
  regional_data <- fread(regional_comparison_path)
  
  # Identify metrics present in both regions
  both_regions <- regional_data %>%
    group_by(Metric, Lag) %>%
    filter(n() == 2) %>%
    arrange(Metric, Lag) %>%
    ungroup()
  
  if (nrow(both_regions) == 0) {
    cat("No metrics significant in both regions. Cannot compute differences.\n")
    return(NULL)
  }
  
  # Calculate effect size differences
  comparison_stats <- both_regions %>%
    group_by(Metric, Lag) %>%
    summarise(
      Bransfield_Coef = Mean_Coefficient[Region == "Bransfield"],
      Central_Coef = Mean_Coefficient[Region == "Central_WAP"],
      Absolute_Diff = abs(Bransfield_Coef - Central_Coef),
      Relative_Diff_Pct = abs(Absolute_Diff / mean(c(Bransfield_Coef, Central_Coef))) * 100,
      Bransfield_N_Models = N_Models[Region == "Bransfield"],
      Central_N_Models = N_Models[Region == "Central_WAP"],
      N_Models_Ratio = max(Bransfield_N_Models, Central_N_Models) / min(Bransfield_N_Models, Central_N_Models),
      Bransfield_N_Obs = Mean_N_Obs[Region == "Bransfield"],
      Central_N_Obs = Mean_N_Obs[Region == "Central_WAP"],
      Sample_Size_Ratio = max(Bransfield_N_Obs, Central_N_Obs) / min(Bransfield_N_Obs, Central_N_Obs),
      Same_Sign = sign(Bransfield_Coef) == sign(Central_Coef),
      .groups = "drop"
    ) %>%
    arrange(desc(Relative_Diff_Pct))
  
  # Summary statistics
  cat("\n=== EFFECT SIZE COMPARISON ===\n\n")
  
  cat(sprintf("Metrics significant in both regions: %d\n", nrow(comparison_stats)))
  cat(sprintf("Metrics with same sign: %d (%.1f%%)\n", 
              sum(comparison_stats$Same_Sign), 
              100 * mean(comparison_stats$Same_Sign)))
  cat(sprintf("\nMean relative difference: %.1f%%\n", mean(comparison_stats$Relative_Diff_Pct)))
  cat(sprintf("Median relative difference: %.1f%%\n", median(comparison_stats$Relative_Diff_Pct)))
  cat(sprintf("Max relative difference: %.1f%%\n", max(comparison_stats$Relative_Diff_Pct)))
  cat(sprintf("Differences < 20%%: %d of %d (%.1f%%)\n",
              sum(comparison_stats$Relative_Diff_Pct < 20),
              nrow(comparison_stats),
              100 * mean(comparison_stats$Relative_Diff_Pct < 20)))
  
  cat("\n=== SAMPLE SIZE BALANCE ===\n\n")
  cat(sprintf("Mean sample size ratio: %.2f:1\n", mean(comparison_stats$Sample_Size_Ratio)))
  cat(sprintf("Max sample size ratio: %.2f:1\n", max(comparison_stats$Sample_Size_Ratio)))
  cat(sprintf("Well-balanced (< 2:1): %d of %d (%.1f%%)\n",
              sum(comparison_stats$Sample_Size_Ratio < 2),
              nrow(comparison_stats),
              100 * mean(comparison_stats$Sample_Size_Ratio < 2)))
  
  cat("\n=== SPATIAL SCALE ROBUSTNESS (N_Models) ===\n\n")
  cat(sprintf("Mean N_Models ratio: %.2f:1\n", mean(comparison_stats$N_Models_Ratio)))
  cat(sprintf("Max N_Models ratio: %.2f:1\n", max(comparison_stats$N_Models_Ratio)))
  cat(sprintf("Similar robustness (< 2:1): %d of %d (%.1f%%)\n",
              sum(comparison_stats$N_Models_Ratio < 2),
              nrow(comparison_stats),
              100 * mean(comparison_stats$N_Models_Ratio < 2)))
  
  # Detailed table
  cat("\n=== DETAILED COMPARISONS ===\n\n")
  
  detailed_output <- comparison_stats %>%
    select(Metric, Lag, Bransfield_Coef, Central_Coef, 
           Absolute_Diff, Relative_Diff_Pct, Same_Sign) %>%
    mutate(
      Bransfield_Coef = round(Bransfield_Coef, 4),
      Central_Coef = round(Central_Coef, 4),
      Absolute_Diff = round(Absolute_Diff, 4),
      Relative_Diff_Pct = round(Relative_Diff_Pct, 1)
    )
  
  print(detailed_output)
  
  # Load interaction results if available
  interaction_summary <- NULL
  if (!is.null(interaction_results_path) && file.exists(interaction_results_path)) {
    cat("\n\n=== INTERACTION TEST RESULTS ===\n\n")
    
    interaction_data <- fread(interaction_results_path)
    
    if (nrow(interaction_data) > 0) {
      interaction_summary <- interaction_data %>%
        mutate(Significant = Bonferroni_Adjusted_pValue < 0.05) %>%
        group_by(Significant) %>%
        summarise(
          N_Tests = n(),
          Mean_P = mean(Interaction_pValue, na.rm = TRUE),
          Median_P = median(Interaction_pValue, na.rm = TRUE),
          .groups = "drop"
        )
      
      cat(sprintf("Total interaction tests: %d\n", nrow(interaction_data)))
      cat(sprintf("Significant interactions (p < 0.05): %d (%.1f%%)\n",
                  sum(interaction_data$Bonferroni_Adjusted_pValue < 0.05, na.rm = TRUE),
                  100 * mean(interaction_data$Bonferroni_Adjusted_pValue < 0.05, na.rm = TRUE)))
      cat(sprintf("Median interaction p-value: %.3f\n", 
                  median(interaction_data$Interaction_pValue, na.rm = TRUE)))
      
      if (sum(interaction_data$Bonferroni_Adjusted_pValue < 0.05, na.rm = TRUE) > 0) {
        cat("\nSignificant interactions found:\n")
        sig_interactions <- interaction_data %>%
          filter(Bonferroni_Adjusted_pValue < 0.05) %>%
          select(Metric, HomeRangeSize, Threshold, Lag, 
                 Interaction_Coefficient, Interaction_pValue, 
                 Bonferroni_Adjusted_pValue) %>%
          arrange(Bonferroni_Adjusted_pValue)
        print(head(sig_interactions, 10))
      }
    }
  }
  
  # Generate recommendation
  cat("\n\n", rep("=", 70), "\n")
  cat("POOLING RECOMMENDATION\n")
  cat(rep("=", 70), "\n\n")
  
  # Criteria for pooling
  criteria_met <- c(
    effect_sizes_similar = mean(comparison_stats$Relative_Diff_Pct < 20) > 0.7,
    same_direction = mean(comparison_stats$Same_Sign) > 0.8,
    balanced_samples = mean(comparison_stats$Sample_Size_Ratio < 2) > 0.7,
    no_significant_interactions = if(!is.null(interaction_summary)) {
      sum(interaction_data$Bonferroni_Adjusted_pValue < 0.05, na.rm = TRUE) == 0
    } else {
      NA
    }
  )
  
  cat("Pooling Criteria:\n")
  cat(sprintf("  ✓ Effect sizes similar (<20%% diff): %s (%.1f%% of comparisons)\n",
              ifelse(criteria_met["effect_sizes_similar"], "YES", "NO"),
              100 * mean(comparison_stats$Relative_Diff_Pct < 20)))
  cat(sprintf("  ✓ Consistent direction (same sign): %s (%.1f%% of comparisons)\n",
              ifelse(criteria_met["same_direction"], "YES", "NO"),
              100 * mean(comparison_stats$Same_Sign)))
  cat(sprintf("  ✓ Balanced sample sizes (<2:1): %s (%.1f%% of comparisons)\n",
              ifelse(criteria_met["balanced_samples"], "YES", "NO"),
              100 * mean(comparison_stats$Sample_Size_Ratio < 2)))
  if (!is.na(criteria_met["no_significant_interactions"])) {
    cat(sprintf("  ✓ No significant interactions: %s\n",
                ifelse(criteria_met["no_significant_interactions"], "YES", "NO")))
  }
  
  recommendation <- if(sum(criteria_met, na.rm = TRUE) >= 3) {
    "POOL REGIONS"
  } else {
    "REPORT SEPARATELY"
  }
  
  cat(sprintf("\n>>> RECOMMENDATION: %s <<<\n", recommendation))
  
  if (recommendation == "POOL REGIONS") {
    cat("\nJustification:\n")
    cat("- Regional effect sizes differ by less than 20% on average\n")
    cat("- Effects show consistent direction across regions\n")
    cat("- Sample sizes are balanced\n")
    if (!is.na(criteria_met["no_significant_interactions"]) && 
        criteria_met["no_significant_interactions"]) {
      cat("- Formal interaction tests found no significant regional differences\n")
    }
    cat("\nPooling will increase statistical power while controlling for\n")
    cat("colony-specific autocorrelation via AR1 structure.\n")
  } else {
    cat("\nCaution: Some criteria suggest regional heterogeneity.\n")
    cat("Consider reporting regions separately with appropriate caveats.\n")
  }
  
  # Save outputs if directory provided
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Save detailed comparison
    write.csv(comparison_stats, 
              file.path(output_dir, "regional_effect_size_comparison.csv"),
              row.names = FALSE)
    
    # Save summary report
    sink(file.path(output_dir, "pooling_justification_report.txt"))
    cat("REGIONAL POOLING JUSTIFICATION REPORT\n")
    cat(rep("=", 70), "\n\n")
    cat("Generated:", Sys.time(), "\n\n")
    
    cat("SUMMARY STATISTICS\n")
    cat(rep("-", 70), "\n")
    cat(sprintf("Metrics in both regions: %d\n", nrow(comparison_stats)))
    cat(sprintf("Mean relative difference: %.1f%%\n", mean(comparison_stats$Relative_Diff_Pct)))
    cat(sprintf("Metrics < 20%% difference: %.1f%%\n", 
                100 * mean(comparison_stats$Relative_Diff_Pct < 20)))
    cat(sprintf("Same direction: %.1f%%\n", 100 * mean(comparison_stats$Same_Sign)))
    cat(sprintf("\nRECOMMENDATION: %s\n", recommendation))
    sink()
    
    cat(sprintf("\nOutputs saved to: %s\n", output_dir))
    cat("  - regional_effect_size_comparison.csv\n")
    cat("  - pooling_justification_report.txt\n")
  }
  
  cat("\n", rep("=", 70), "\n\n")
  
  # Return results
  return(list(
    comparison_stats = comparison_stats,
    interaction_summary = interaction_summary,
    criteria_met = criteria_met,
    recommendation = recommendation
  ))
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if (FALSE) {  # Set to TRUE to run
  
  results_path <- "C:/Users/michael.wethington.BRILOON/OneDrive - Biodiversity Research Institute/Documents/Manuscripts - Antarctica/FrostBound_AQ/RStudioProject/pipelines/gentoo_abundance_analysis/results"
  
  pooling_analysis <- analyze_pooling_justification(
    regional_comparison_path = file.path(results_path, "Regional_Comparison_Summary.csv"),
    interaction_results_path = file.path(results_path, "Interaction_Analysis", "interaction_tests_all.csv"),
    output_dir = file.path(results_path, "Pooling_Justification")
  )
  
  # Access results
  print(pooling_analysis$recommendation)
  print(pooling_analysis$comparison_stats)
}