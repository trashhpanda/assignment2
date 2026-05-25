############################### Task 2: Bootstrap on regression model ###############################

# Create output directories if they don't exist
dir.create("results", showWarnings = FALSE)
dir.create("figure", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(123)

# Load libraries
library(boot)
library(MASS)

# Load Animals dataset
data(Animals)

# Filter: keep only observations where body > 0 AND brain > 0 and no missing values
animals_clean <- Animals[
  Animals$body > 0 & Animals$brain > 0 & complete.cases(Animals[, c("body", "brain")]),
]

# Create log-transformed variables
animals_clean$log_body <- log(animals_clean$body)
animals_clean$log_brain <- log(animals_clean$brain)

# Fit the model on the full sample
model_original <- lm(log_body ~ log_brain, data = animals_clean)

# Extract R² from original sample
r2_original <- summary(model_original)$r.squared

print(paste("Original R²:", round(r2_original, 4)))

############################### Bootstrap statistic function ###############################

# Bootstrap statistic function for boot package
# data: the dataset
# indices: bootstrap indices provided by boot()
# Returns: R² value
r2_statistic <- function(data, indices) {
  # Select bootstrap sample using provided indices
  boot_sample <- data[indices, ]
  
  # Fit regression model on bootstrap sample
  boot_model <- lm(log_body ~ log_brain, data = boot_sample)
  
  # Calculate in-sample R²
  boot_r2 <- summary(boot_model)$r.squared
  
  return(boot_r2)
}

############################### Run bootstrap ###############################

# Run bootstrap with 10,000 iterations
bootstrap_results <- boot(data = animals_clean, 
                         statistic = r2_statistic, 
                         R = 10000)

# Extract bootstrap R² estimates
r2_bootstrap_estimates <- bootstrap_results$t[, 1]

# Calculate mean of bootstrap estimates
r2_bootstrap_mean <- mean(r2_bootstrap_estimates)

# Calculate bootstrap bias
bootstrap_bias <- r2_bootstrap_mean - r2_original

# Calculate 95% confidence interval using percentile method
ci_percentile <- boot.ci(bootstrap_results, type = "perc", conf = 0.95)

# Extract CI bounds
ci_lower <- ci_percentile$percent[4]
ci_upper <- ci_percentile$percent[5]

############################### Output and reporting ###############################

cat("\n=== Task 2: Bootstrap on Regression Model ===\n\n")

cat("Data Summary:\n")
cat("- Sample size (after filtering):", nrow(animals_clean), "\n")
cat("- Model: log(body) ~ log(brain)\n")
cat("- Bootstrap type: nonparametric case bootstrap\n")
cat("- Number of bootstrap iterations: 10,000\n\n")

cat("Results:\n")
cat("- Original R² (R²_original):", round(r2_original, 4), "\n")
cat("- Bootstrap mean R² (R²_bootstrap_mean):", round(r2_bootstrap_mean, 4), "\n")
cat("- Bootstrap bias:", round(bootstrap_bias, 6), "\n")
cat("- 95% Percentile CI for R²: [", round(ci_lower, 4), ", ", round(ci_upper, 4), "]\n\n")

cat("Interpretation:\n")
# Check if original R² lies inside the CI
inside_ci <- (r2_original >= ci_lower) & (r2_original <= ci_upper)
if (inside_ci) {
  cat("- The original R² lies INSIDE the bootstrap confidence interval.\n")
} else {
  cat("- The original R² lies OUTSIDE the bootstrap confidence interval.\n")
}

# Comment on bias magnitude
if (abs(bootstrap_bias) < 0.01) {
  cat("- The bootstrap bias is small (", round(bootstrap_bias, 6), 
      ") relative to the original R² estimate.\n")
} else {
  cat("- The bootstrap bias is non-negligible (", round(bootstrap_bias, 6), 
      ") relative to the original R² estimate.\n")
}

cat("- The confidence interval provides a range of plausible R² values based on\n")
cat("  resampling the observed data, reflecting uncertainty in model fit.\n\n")

############################### Save results ###############################

task2_results <- list(
  r2_original = r2_original,
  r2_bootstrap_mean = r2_bootstrap_mean,
  bootstrap_bias = bootstrap_bias,
  ci_lower = ci_lower,
  ci_upper = ci_upper,
  r2_estimates = r2_bootstrap_estimates,
  bootstrap_object = bootstrap_results,
  sample_size = nrow(animals_clean)
)

saveRDS(task2_results, file = "results/task2_results.rds")
cat("Results saved to results/task2_results.rds\n\n")

############################### Histogram and visualization ###############################

# Create histogram of bootstrap R² estimates
png(filename = "figure/task2_histogram.png", width = 800, height = 600)

hist(r2_bootstrap_estimates, 
     main = "Bootstrap R² Distribution (N = 10,000 iterations)",
     xlab = "R²",
     ylab = "Frequency",
     col = "lightblue",
     border = "darkblue",
     breaks = 50)

# Add vertical lines
abline(v = r2_original, col = "red", lwd = 2, lty = 2)
abline(v = r2_bootstrap_mean, col = "green", lwd = 2, lty = 2)
abline(v = ci_lower, col = "purple", lwd = 2, lty = 3)
abline(v = ci_upper, col = "purple", lwd = 2, lty = 3)

legend("topright", 
       legend = c("R²_original", "R²_bootstrap_mean", "95% CI bounds"),
       col = c("red", "green", "purple"),
       lty = c(2, 2, 3),
       lwd = 2)

dev.off()

cat("Histogram saved to figure/task2_histogram.png\n")
