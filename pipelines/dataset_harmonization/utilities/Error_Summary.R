# Sort the heatmap data by SMAPE in descending order
top_errors <- heatmap_data %>%
  arrange(desc(SMAPE)) %>%
  head(10)  # Modify '10' to the number of top bins you want to examine

# View the top bins contributing the most error
print(top_errors)

# Calculate total or average SMAPE for specific binned ranges if needed:
# For example, summarizing error for AMSR values between 0.0 and 0.2
summary_error <- heatmap_data %>%
  filter(AMSR_bin %in% c("(0,0.05]", "(0.05,0.1]", "(0.1,0.15]", "(0.15,0.2]")) %>%
  summarize(Total_Error = sum(SMAPE), Avg_Error = mean(SMAPE))

print(summary_error)
