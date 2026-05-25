############################### Task 3: Bootstrap salary data analysis ###############################

# Create output directories if they don't exist
dir.create("results", showWarnings = FALSE)
dir.create("figure", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(123)

# Load libraries
library(boot)
library(wooldridge)

# Load nbasal dataset
data(nbasal)

# Display column names to verify structure
cat("Column names in nbasal:\n")
print(names(nbasal))
cat("\n")

# Note: In wooldridge::nbasal, marital status is stored in variable 'marr' (not 'married')
# wage is measured in thousands of USD

# Remove observations with missing values in wage or marr
nbasal_clean <- nbasal[complete.cases(nbasal[, c("wage", "marr")]), ]

cat("\n=== Task 3: Bootstrap Salary Data Analysis ===\n")
cat("(Wage is measured in thousands of USD)\n\n")

cat("Data Summary:\n")
cat("- Total observations (after removing NA):", nrow(nbasal_clean), "\n")
cat("- Bootstrap iterations: 10,000\n")
cat("- Bootstrap type: nonparametric case bootstrap\n")
cat("- Confidence level: 95%\n\n")

############################### PART 1: Whole sample analysis ###############################

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("PART 1: WHOLE SAMPLE ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Define bootstrap statistics for mean and median (whole sample)
mean_statistic <- function(data, indices) {
  return(mean(data[indices]))
}

median_statistic <- function(data, indices) {
  return(median(data[indices]))
}

# Bootstrap for mean (whole sample)
boot_mean_whole <- boot(data = nbasal_clean$wage, 
                        statistic = mean_statistic, 
                        R = 10000)

# Bootstrap for median (whole sample)
boot_median_whole <- boot(data = nbasal_clean$wage, 
                          statistic = median_statistic, 
                          R = 10000)

# Original estimates
mean_original <- mean(nbasal_clean$wage)
median_original <- median(nbasal_clean$wage)

# Bootstrap means
mean_bootstrap_mean <- mean(boot_mean_whole$t[, 1])
median_bootstrap_mean <- mean(boot_median_whole$t[, 1])

# Bootstrap bias
mean_bias <- mean_bootstrap_mean - mean_original
median_bias <- median_bootstrap_mean - median_original

# Confidence intervals using percentile method
ci_mean_whole <- boot.ci(boot_mean_whole, type = "perc", conf = 0.95)
ci_median_whole <- boot.ci(boot_median_whole, type = "perc", conf = 0.95)

mean_ci_lower <- ci_mean_whole$percent[4]
mean_ci_upper <- ci_mean_whole$percent[5]
median_ci_lower <- ci_median_whole$percent[4]
median_ci_upper <- ci_median_whole$percent[5]

cat("Whole sample estimates:\n")
cat("- Original mean wage:", round(mean_original, 2), "\n")
cat("- Bootstrap mean estimate:", round(mean_bootstrap_mean, 2), "\n")
cat("- Bootstrap bias:", round(mean_bias, 4), "\n")
cat("- 95% CI for mean: [", round(mean_ci_lower, 2), ", ", round(mean_ci_upper, 2), "]\n\n")

cat("- Original median wage:", round(median_original, 2), "\n")
cat("- Bootstrap median estimate:", round(median_bootstrap_mean, 2), "\n")
cat("- Bootstrap bias:", round(median_bias, 4), "\n")
cat("- 95% CI for median: [", round(median_ci_lower, 2), ", ", round(median_ci_upper, 2), "]\n\n")

############################### PART 2: Married vs unmarried ###############################

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("PART 2: MARRIED VS UNMARRIED COMPARISON\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Split data by marital status
married_data <- nbasal_clean[nbasal_clean$marr == 1, ]
unmarried_data <- nbasal_clean[nbasal_clean$marr == 0, ]

cat("Group summary:\n")
cat("- Married (married = 1): n =", nrow(married_data), "\n")
cat("- Unmarried (married = 0): n =", nrow(unmarried_data), "\n\n")

# Bootstrap for each group
boot_mean_married <- boot(data = married_data$wage, 
                          statistic = mean_statistic, 
                          R = 10000)
boot_median_married <- boot(data = married_data$wage, 
                            statistic = median_statistic, 
                            R = 10000)

boot_mean_unmarried <- boot(data = unmarried_data$wage, 
                             statistic = mean_statistic, 
                             R = 10000)
boot_median_unmarried <- boot(data = unmarried_data$wage, 
                               statistic = median_statistic, 
                               R = 10000)

# Married group estimates
mean_original_married <- mean(married_data$wage)
median_original_married <- median(married_data$wage)
mean_bootstrap_married <- mean(boot_mean_married$t[, 1])
median_bootstrap_married <- mean(boot_median_married$t[, 1])
mean_bias_married <- mean_bootstrap_married - mean_original_married
median_bias_married <- median_bootstrap_married - median_original_married

ci_mean_married <- boot.ci(boot_mean_married, type = "perc", conf = 0.95)
ci_median_married <- boot.ci(boot_median_married, type = "perc", conf = 0.95)

# Unmarried group estimates
mean_original_unmarried <- mean(unmarried_data$wage)
median_original_unmarried <- median(unmarried_data$wage)
mean_bootstrap_unmarried <- mean(boot_mean_unmarried$t[, 1])
median_bootstrap_unmarried <- mean(boot_median_unmarried$t[, 1])
mean_bias_unmarried <- mean_bootstrap_unmarried - mean_original_unmarried
median_bias_unmarried <- median_bootstrap_unmarried - median_original_unmarried

ci_mean_unmarried <- boot.ci(boot_mean_unmarried, type = "perc", conf = 0.95)
ci_median_unmarried <- boot.ci(boot_median_unmarried, type = "perc", conf = 0.95)

cat("MARRIED (married = 1):\n")
cat("- Original mean wage:", round(mean_original_married, 2), "\n")
cat("- Bootstrap mean estimate:", round(mean_bootstrap_married, 2), "\n")
cat("- Bootstrap bias:", round(mean_bias_married, 4), "\n")
cat("- 95% CI for mean: [", round(ci_mean_married$percent[4], 2), ", ", round(ci_mean_married$percent[5], 2), "]\n\n")

cat("- Original median wage:", round(median_original_married, 2), "\n")
cat("- Bootstrap median estimate:", round(median_bootstrap_married, 2), "\n")
cat("- Bootstrap bias:", round(median_bias_married, 4), "\n")
cat("- 95% CI for median: [", round(ci_median_married$percent[4], 2), ", ", round(ci_median_married$percent[5], 2), "]\n\n")

cat("UNMARRIED (married = 0):\n")
cat("- Original mean wage:", round(mean_original_unmarried, 2), "\n")
cat("- Bootstrap mean estimate:", round(mean_bootstrap_unmarried, 2), "\n")
cat("- Bootstrap bias:", round(mean_bias_unmarried, 4), "\n")
cat("- 95% CI for mean: [", round(ci_mean_unmarried$percent[4], 2), ", ", round(ci_mean_unmarried$percent[5], 2), "]\n\n")

cat("- Original median wage:", round(median_original_unmarried, 2), "\n")
cat("- Bootstrap median estimate:", round(median_bootstrap_unmarried, 2), "\n")
cat("- Bootstrap bias:", round(median_bias_unmarried, 4), "\n")
cat("- 95% CI for median: [", round(ci_median_unmarried$percent[4], 2), ", ", round(ci_median_unmarried$percent[5], 2), "]\n\n")

# Descriptive comparison of confidence intervals
cat("Overlapping confidence intervals (descriptive analysis):\n")
if (ci_mean_married$percent[4] > ci_mean_unmarried$percent[5] | ci_mean_unmarried$percent[4] > ci_mean_married$percent[5]) {
  cat("- Mean CI for married and unmarried DO NOT overlap.\n")
} else {
  cat("- Mean CI for married and unmarried DO overlap.\n")
}

if (ci_median_married$percent[4] > ci_median_unmarried$percent[5] | ci_median_unmarried$percent[4] > ci_median_married$percent[5]) {
  cat("- Median CI for married and unmarried DO NOT overlap.\n")
} else {
  cat("- Median CI for married and unmarried DO overlap.\n")
}

cat("\nNote: Overlapping confidence intervals alone are not a formal test of equality.\n")
cat("Therefore, we additionally analyse the bootstrap confidence interval for the difference.\n\n")

############################### PART 3: Differences between groups ###############################

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("PART 3: DIFFERENCE IN MEANS AND MEDIANS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Bootstrap function for difference in means
diff_mean_statistic <- function(data, indices) {
  boot_data <- data[indices, ]
  
  mean_married <- mean(boot_data$wage[boot_data$marr == 1])
  mean_unmarried <- mean(boot_data$wage[boot_data$marr == 0])
  
  return(mean_married - mean_unmarried)
}

# Bootstrap function for difference in medians
diff_median_statistic <- function(data, indices) {
  boot_data <- data[indices, ]
  
  median_married <- median(boot_data$wage[boot_data$marr == 1])
  median_unmarried <- median(boot_data$wage[boot_data$marr == 0])
  
  return(median_married - median_unmarried)
}

# Bootstrap differences using stratified bootstrap via strata argument
boot_diff_mean <- boot(
  data = nbasal_clean,
  statistic = diff_mean_statistic,
  R = 10000,
  strata = nbasal_clean$marr
)

boot_diff_median <- boot(
  data = nbasal_clean,
  statistic = diff_median_statistic,
  R = 10000,
  strata = nbasal_clean$marr
)

# Original differences
diff_mean_original <- mean_original_married - mean_original_unmarried
diff_median_original <- median_original_married - median_original_unmarried

# Bootstrap estimates of differences
diff_mean_bootstrap <- mean(boot_diff_mean$t[, 1])
diff_median_bootstrap <- mean(boot_diff_median$t[, 1])

# Bias of differences
diff_mean_bias <- diff_mean_bootstrap - diff_mean_original
diff_median_bias <- diff_median_bootstrap - diff_median_original

# Confidence intervals for differences
ci_diff_mean <- boot.ci(boot_diff_mean, type = "perc", conf = 0.95)
ci_diff_median <- boot.ci(boot_diff_median, type = "perc", conf = 0.95)

diff_mean_ci_lower <- ci_diff_mean$percent[4]
diff_mean_ci_upper <- ci_diff_mean$percent[5]
diff_median_ci_lower <- ci_diff_median$percent[4]
diff_median_ci_upper <- ci_diff_median$percent[5]

cat("Difference in means (married - unmarried):\n")
cat("- Original difference:", round(diff_mean_original, 2), "\n")
cat("- Bootstrap mean difference:", round(diff_mean_bootstrap, 2), "\n")
cat("- Bootstrap bias:", round(diff_mean_bias, 4), "\n")
cat("- 95% CI for difference: [", round(diff_mean_ci_lower, 2), ", ", round(diff_mean_ci_upper, 2), "]\n")

# Check if CI contains 0
if (diff_mean_ci_lower <= 0 & diff_mean_ci_upper >= 0) {
  cat("- Interpretation: CI contains 0. Bootstrap does not provide evidence of a difference in means.\n\n")
  conclusion_mean <- "No evidence of difference"
} else {
  cat("- Interpretation: CI does NOT contain 0. Bootstrap provides evidence of a difference in means.\n\n")
  conclusion_mean <- "Evidence of difference"
}

cat("Difference in medians (married - unmarried):\n")
cat("- Original difference:", round(diff_median_original, 2), "\n")
cat("- Bootstrap mean difference:", round(diff_median_bootstrap, 2), "\n")
cat("- Bootstrap bias:", round(diff_median_bias, 4), "\n")
cat("- 95% CI for difference: [", round(diff_median_ci_lower, 2), ", ", round(diff_median_ci_upper, 2), "]\n")

# Check if CI contains 0
if (diff_median_ci_lower <= 0 & diff_median_ci_upper >= 0) {
  cat("- Interpretation: CI contains 0. Bootstrap does not provide evidence of a difference in medians.\n\n")
  conclusion_median <- "No evidence of difference"
} else {
  cat("- Interpretation: CI does NOT contain 0. Bootstrap provides evidence of a difference in medians.\n\n")
  conclusion_median <- "Evidence of difference"
}

############################### PART 4: Optional - Shortest confidence interval ###############################

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("PART 4: SHORTEST CONFIDENCE INTERVAL (Optional)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Function to compute shortest confidence interval
shortest_ci <- function(boot_estimates, conf = 0.95) {
  n_boot <- length(boot_estimates)
  window_size <- ceiling(n_boot * conf)
  
  sorted_estimates <- sort(boot_estimates)
  
  # Find all possible windows
  n_windows <- n_boot - window_size + 1
  widths <- numeric(n_windows)
  
  for (i in 1:n_windows) {
    lower <- sorted_estimates[i]
    upper <- sorted_estimates[i + window_size - 1]
    widths[i] <- upper - lower
  }
  
  # Find window with minimum width
  min_idx <- which.min(widths)
  shortest_lower <- sorted_estimates[min_idx]
  shortest_upper <- sorted_estimates[min_idx + window_size - 1]
  shortest_width <- shortest_upper - shortest_lower
  
  return(list(lower = shortest_lower, upper = shortest_upper, width = shortest_width))
}

# Compute shortest CI for mean and median
shortest_mean_whole <- shortest_ci(boot_mean_whole$t[, 1], conf = 0.95)
shortest_median_whole <- shortest_ci(boot_median_whole$t[, 1], conf = 0.95)

# Compare with percentile CI
percentile_mean_width <- mean_ci_upper - mean_ci_lower
percentile_median_width <- median_ci_upper - median_ci_lower

cat("WHOLE SAMPLE - MEAN:\n")
cat("- Percentile CI: [", round(mean_ci_lower, 2), ", ", round(mean_ci_upper, 2), "], width =", round(percentile_mean_width, 2), "\n")
cat("- Shortest CI:  [", round(shortest_mean_whole$lower, 2), ", ", round(shortest_mean_whole$upper, 2), "], width =", round(shortest_mean_whole$width, 2), "\n\n")

cat("WHOLE SAMPLE - MEDIAN:\n")
cat("- Percentile CI: [", round(median_ci_lower, 2), ", ", round(median_ci_upper, 2), "], width =", round(percentile_median_width, 2), "\n")
cat("- Shortest CI:  [", round(shortest_median_whole$lower, 2), ", ", round(shortest_median_whole$upper, 2), "], width =", round(shortest_median_whole$width, 2), "\n\n")

############################### Save results ###############################

task3_results <- list(
  # Whole sample
  mean_original = mean_original,
  mean_bootstrap_mean = mean_bootstrap_mean,
  mean_bias = mean_bias,
  mean_ci_lower = mean_ci_lower,
  mean_ci_upper = mean_ci_upper,
  median_original = median_original,
  median_bootstrap_mean = median_bootstrap_mean,
  median_bias = median_bias,
  median_ci_lower = median_ci_lower,
  median_ci_upper = median_ci_upper,
  
  # Married group
  n_married = nrow(married_data),
  mean_original_married = mean_original_married,
  mean_bootstrap_married = mean_bootstrap_married,
  mean_bias_married = mean_bias_married,
  mean_ci_married_lower = ci_mean_married$percent[4],
  mean_ci_married_upper = ci_mean_married$percent[5],
  median_original_married = median_original_married,
  median_bootstrap_married = median_bootstrap_married,
  median_bias_married = median_bias_married,
  median_ci_married_lower = ci_median_married$percent[4],
  median_ci_married_upper = ci_median_married$percent[5],
  
  # Unmarried group
  n_unmarried = nrow(unmarried_data),
  mean_original_unmarried = mean_original_unmarried,
  mean_bootstrap_unmarried = mean_bootstrap_unmarried,
  mean_bias_unmarried = mean_bias_unmarried,
  mean_ci_unmarried_lower = ci_mean_unmarried$percent[4],
  mean_ci_unmarried_upper = ci_mean_unmarried$percent[5],
  median_original_unmarried = median_original_unmarried,
  median_bootstrap_unmarried = median_bootstrap_unmarried,
  median_bias_unmarried = median_bias_unmarried,
  median_ci_unmarried_lower = ci_median_unmarried$percent[4],
  median_ci_unmarried_upper = ci_median_unmarried$percent[5],
  
  # Differences
  diff_mean_original = diff_mean_original,
  diff_mean_bootstrap = diff_mean_bootstrap,
  diff_mean_bias = diff_mean_bias,
  diff_mean_ci_lower = diff_mean_ci_lower,
  diff_mean_ci_upper = diff_mean_ci_upper,
  diff_mean_conclusion = conclusion_mean,
  diff_median_original = diff_median_original,
  diff_median_bootstrap = diff_median_bootstrap,
  diff_median_bias = diff_median_bias,
  diff_median_ci_lower = diff_median_ci_lower,
  diff_median_ci_upper = diff_median_ci_upper,
  diff_median_conclusion = conclusion_median,
  
  # Bootstrap objects for visualization
  boot_mean_whole = boot_mean_whole,
  boot_median_whole = boot_median_whole,
  boot_diff_mean = boot_diff_mean,
  boot_diff_median = boot_diff_median
)

saveRDS(task3_results, file = "results/task3_results.rds")
cat("Results saved to results/task3_results.rds\n\n")

############################### Histograms and visualizations ###############################

# Set up for multiple plots
png(filename = "figure/task3_whole_sample_mean.png", width = 800, height = 600)
hist(boot_mean_whole$t[, 1],
     main = "Bootstrap Distribution of Mean Wage (Whole Sample)",
     xlab = "Mean Wage",
     ylab = "Frequency",
     col = "lightblue",
     border = "darkblue",
     breaks = 50)
abline(v = mean_original, col = "red", lwd = 2, lty = 2)
abline(v = mean_bootstrap_mean, col = "green", lwd = 2, lty = 2)
abline(v = mean_ci_lower, col = "purple", lwd = 2, lty = 3)
abline(v = mean_ci_upper, col = "purple", lwd = 2, lty = 3)
legend("topright",
       legend = c("Original", "Bootstrap mean", "95% CI bounds"),
       col = c("red", "green", "purple"),
       lty = c(2, 2, 3),
       lwd = 2)
dev.off()
cat("Saved: figure/task3_whole_sample_mean.png\n")

png(filename = "figure/task3_whole_sample_median.png", width = 800, height = 600)
hist(boot_median_whole$t[, 1],
     main = "Bootstrap Distribution of Median Wage (Whole Sample)",
     xlab = "Median Wage",
     ylab = "Frequency",
     col = "lightblue",
     border = "darkblue",
     breaks = 50)
abline(v = median_original, col = "red", lwd = 2, lty = 2)
abline(v = median_bootstrap_mean, col = "green", lwd = 2, lty = 2)
abline(v = median_ci_lower, col = "purple", lwd = 2, lty = 3)
abline(v = median_ci_upper, col = "purple", lwd = 2, lty = 3)
legend("topright",
       legend = c("Original", "Bootstrap mean", "95% CI bounds"),
       col = c("red", "green", "purple"),
       lty = c(2, 2, 3),
       lwd = 2)
dev.off()
cat("Saved: figure/task3_whole_sample_median.png\n")

png(filename = "figure/task3_married_mean.png", width = 800, height = 600)
hist(boot_mean_married$t[, 1],
     main = "Bootstrap Distribution of Mean Wage (Married = 1)",
     xlab = "Mean Wage",
     ylab = "Frequency",
     col = "lightcyan",
     border = "darkcyan",
     breaks = 50)
abline(v = mean_original_married, col = "red", lwd = 2, lty = 2)
abline(v = mean_bootstrap_married, col = "green", lwd = 2, lty = 2)
abline(v = ci_mean_married$percent[4], col = "purple", lwd = 2, lty = 3)
abline(v = ci_mean_married$percent[5], col = "purple", lwd = 2, lty = 3)
legend("topright",
       legend = c("Original", "Bootstrap mean", "95% CI bounds"),
       col = c("red", "green", "purple"),
       lty = c(2, 2, 3),
       lwd = 2)
dev.off()
cat("Saved: figure/task3_married_mean.png\n")

png(filename = "figure/task3_unmarried_mean.png", width = 800, height = 600)
hist(boot_mean_unmarried$t[, 1],
     main = "Bootstrap Distribution of Mean Wage (Married = 0)",
     xlab = "Mean Wage",
     ylab = "Frequency",
     col = "lightgoldenrodyellow",
     border = "goldenrod",
     breaks = 50)
abline(v = mean_original_unmarried, col = "red", lwd = 2, lty = 2)
abline(v = mean_bootstrap_unmarried, col = "green", lwd = 2, lty = 2)
abline(v = ci_mean_unmarried$percent[4], col = "purple", lwd = 2, lty = 3)
abline(v = ci_mean_unmarried$percent[5], col = "purple", lwd = 2, lty = 3)
legend("topright",
       legend = c("Original", "Bootstrap mean", "95% CI bounds"),
       col = c("red", "green", "purple"),
       lty = c(2, 2, 3),
       lwd = 2)
dev.off()
cat("Saved: figure/task3_unmarried_mean.png\n")

png(filename = "figure/task3_diff_mean.png", width = 800, height = 600)
hist(boot_diff_mean$t[, 1],
     main = "Bootstrap Distribution of Mean Difference\n(Married - Unmarried)",
     xlab = "Mean Difference",
     ylab = "Frequency",
     col = "lightcoral",
     border = "darkred",
     breaks = 50)
abline(v = 0, col = "black", lwd = 1, lty = 1)
abline(v = diff_mean_original, col = "red", lwd = 2, lty = 2)
abline(v = diff_mean_bootstrap, col = "green", lwd = 2, lty = 2)
abline(v = diff_mean_ci_lower, col = "purple", lwd = 2, lty = 3)
abline(v = diff_mean_ci_upper, col = "purple", lwd = 2, lty = 3)
legend("topright",
       legend = c("Zero", "Original", "Bootstrap mean", "95% CI bounds"),
       col = c("black", "red", "green", "purple"),
       lty = c(1, 2, 2, 3),
       lwd = 2)
dev.off()
cat("Saved: figure/task3_diff_mean.png\n")

png(filename = "figure/task3_diff_median.png", width = 800, height = 600)
hist(boot_diff_median$t[, 1],
     main = "Bootstrap Distribution of Median Difference\n(Married - Unmarried)",
     xlab = "Median Difference",
     ylab = "Frequency",
     col = "lightcoral",
     border = "darkred",
     breaks = 50)
abline(v = 0, col = "black", lwd = 1, lty = 1)
abline(v = diff_median_original, col = "red", lwd = 2, lty = 2)
abline(v = diff_median_bootstrap, col = "green", lwd = 2, lty = 2)
abline(v = diff_median_ci_lower, col = "purple", lwd = 2, lty = 3)
abline(v = diff_median_ci_upper, col = "purple", lwd = 2, lty = 3)
legend("topright",
       legend = c("Zero", "Original", "Bootstrap mean", "95% CI bounds"),
       col = c("black", "red", "green", "purple"),
       lty = c(1, 2, 2, 3),
       lwd = 2)
dev.off()
cat("Saved: figure/task3_diff_median.png\n")

############################### Final Summary ###############################

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("FINAL SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Key findings:\n\n")
cat("1. Mean Salary Comparison:\n")
cat("   - Married players: ", round(mean_original_married, 2), " thousand USD\n", sep = "")
cat("   - Unmarried players: ", round(mean_original_unmarried, 2), " thousand USD\n", sep = "")
cat("   - Difference: ", round(diff_mean_original, 2), " thousand USD\n", sep = "")
cat("   - Bootstrap 95% CI for difference: [", round(diff_mean_ci_lower, 2), ", ", round(diff_mean_ci_upper, 2), "]\n", sep = "")
cat("   - Conclusion: ", conclusion_mean, " in mean salary between groups.\n\n", sep = "")

cat("2. Median Salary Comparison:\n")
cat("   - Married players: ", round(median_original_married, 2), " thousand USD\n", sep = "")
cat("   - Unmarried players: ", round(median_original_unmarried, 2), " thousand USD\n", sep = "")
cat("   - Difference: ", round(diff_median_original, 2), " thousand USD\n", sep = "")
cat("   - Bootstrap 95% CI for difference: [", round(diff_median_ci_lower, 2), ", ", round(diff_median_ci_upper, 2), "]\n", sep = "")
cat("   - Conclusion: ", conclusion_median, " in median salary between groups.\n\n", sep = "")

cat("3. Shortest Confidence Intervals:\n")
cat("   - Mean (whole sample): percentile width = ", round(percentile_mean_width, 2), 
    ", shortest CI width = ", round(shortest_mean_whole$width, 2), "\n", sep = "")
cat("   - Median (whole sample): percentile width = ", round(percentile_median_width, 2), 
    ", shortest CI width = ", round(shortest_median_whole$width, 2), "\n\n", sep = "")

cat("Note: Bootstrap for group differences was stratified by marital status to preserve group structure.\n")
cat("All estimates based on 10,000 bootstrap resamples with 95% confidence level.\n")

cat("\n=== Task 3 Complete ===\n")
