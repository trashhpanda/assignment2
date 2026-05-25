# ============================================================================
# Task 4: Bootstrap Confidence Areas for L3 Lorenz Curves
# ============================================================================

library(boot)
library(wooldridge)

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("figure", showWarnings = FALSE)

# ============================================================================
# PART 1: Gamma(α=2, β=4) Distribution
# ============================================================================

cat("\n=== PART 1: Gamma(α=2, β=4) Distribution ===\n\n")

set.seed(42)

# Generate sample
alpha_true <- 2
beta_true <- 4
n_gamma <- 60
gamma_sample <- rgamma(n = n_gamma, shape = alpha_true, scale = beta_true)

cat("Gamma sample generated: n =", n_gamma, ", α =", alpha_true, ", β =", beta_true, "\n")
cat("Sample mean:", mean(gamma_sample), "\n")
cat("Theoretical mean:", alpha_true * beta_true, "\n\n")

# L3 function - compute Lorenz curve L3(p; Q) for given dataset and p values
compute_l3 <- function(x, p_values) {
  # Q(p/2) and Q(1-p/2) - empirical quantile function (type=1: inverse of empirical CDF)
  Q_lower <- quantile(x, probs = p_values / 2, type = 1, names = FALSE)
  Q_upper <- quantile(x, probs = 1 - p_values / 2, type = 1, names = FALSE)
  
  # L3(p; Q) = p² · Q(p/2) / (Q(p/2) + Q(1-p/2))
  L3 <- (p_values^2 * Q_lower) / (Q_lower + Q_upper)
  
  return(L3)
}

# True L3 function - compute theoretical Lorenz curve using qgamma
compute_l3_gamma_true <- function(p_values, alpha, beta) {
  # Theoretical quantile function from Gamma distribution
  Q_lower <- qgamma(p_values / 2, shape = alpha, scale = beta)
  Q_upper <- qgamma(1 - p_values / 2, shape = alpha, scale = beta)
  
  # L3(p; Q) = p² · Q(p/2) / (Q(p/2) + Q(1-p/2))
  L3_true <- (p_values^2 * Q_lower) / (Q_lower + Q_upper)
  
  return(L3_true)
}

# Bootstrap statistic function - returns L3 curve for bootstrap sample
l3_bootstrap_stat <- function(data, indices, p_values) {
  boot_sample <- data[indices]
  return(compute_l3(boot_sample, p_values))
}

# Evaluation points
p_values_gamma <- seq(0.05, 0.95, length.out = 100)
cat("Number of p values:", length(p_values_gamma), "\n")
cat("p range: [", min(p_values_gamma), ",", max(p_values_gamma), "]\n\n")

# Original L3 curve
L3_original_gamma <- compute_l3(gamma_sample, p_values_gamma)

# True L3 curve
L3_true_gamma <- compute_l3_gamma_true(p_values_gamma, alpha = alpha_true, beta = beta_true)

# Bootstrap L3 curves - boot object returns matrix (B × m)
cat("Running bootstrap (R=10000)...\n")
boot_l3_gamma <- boot(
  data = gamma_sample,
  statistic = l3_bootstrap_stat,
  R = 10000,
  p_values = p_values_gamma
)

cat("Bootstrap complete!\n")
cat("Boot object dimensions: ", dim(boot_l3_gamma$t), " (iterations × p values)\n\n")

# Extract bootstrap statistics for each p
# boot_l3_gamma$t is (10000 × 100) matrix
L3_mean_gamma <- apply(boot_l3_gamma$t, 2, mean)  # mean across bootstrap samples
L3_ci_lower_gamma <- apply(boot_l3_gamma$t, 2, quantile, probs = 0.025)  # 2.5% quantile
L3_ci_upper_gamma <- apply(boot_l3_gamma$t, 2, quantile, probs = 0.975)  # 97.5% quantile

cat("L3 original - range: [", min(L3_original_gamma), ", ", max(L3_original_gamma), "]\n")
cat("L3 true     - range: [", min(L3_true_gamma), ", ", max(L3_true_gamma), "]\n")
cat("L3 mean     - range: [", min(L3_mean_gamma), ", ", max(L3_mean_gamma), "]\n")
cat("L3 CI lower - range: [", min(L3_ci_lower_gamma), ", ", max(L3_ci_lower_gamma), "]\n")
cat("L3 CI upper - range: [", min(L3_ci_upper_gamma), ", ", max(L3_ci_upper_gamma), "]\n\n")

# Check if true L3 curve is inside confidence area
true_inside_ci <- (L3_true_gamma >= L3_ci_lower_gamma) &
                  (L3_true_gamma <= L3_ci_upper_gamma)
inside_count <- sum(true_inside_ci)
inside_pct <- 100 * inside_count / length(p_values_gamma)

cat("True L3 curve inside confidence area at",
    inside_count, "out of", length(p_values_gamma),
    "p values (", round(inside_pct, 2), "%)\n\n")

# Visualization - Part 1
png(file = "figure/task4_part1_gamma.png", width = 8, height = 6, units = "in", res = 150)

# Empty plot first to avoid layering issues
plot(p_values_gamma, L3_original_gamma,
     type = "n",
     main = "L3 Lorenz Curve - Gamma(α=2, β=4) Distribution\nwith 95% Bootstrap Confidence Area",
     xlab = "p", ylab = "L3(p; Q)",
     ylim = c(min(L3_ci_lower_gamma) - 0.05, max(L3_ci_upper_gamma) + 0.05),
     cex.main = 1.2, cex.lab = 1.1)

# Confidence area - polygon
polygon(
  c(p_values_gamma, rev(p_values_gamma)),
  c(L3_ci_lower_gamma, rev(L3_ci_upper_gamma)),
  col = rgb(0.3, 0.6, 1, alpha = 0.25),  # transparent blue
  border = NA
)

# Original empirical L3 curve
lines(p_values_gamma, L3_original_gamma, lwd = 2.5, col = "red", lty = 1)

# True (theoretical) L3 curve
lines(p_values_gamma, L3_true_gamma, lwd = 2.5, col = "purple", lty = 1)

# Bootstrap mean curve
lines(p_values_gamma, L3_mean_gamma, lwd = 2, col = "green", lty = 2)

# CI bounds
lines(p_values_gamma, L3_ci_lower_gamma, lwd = 1.2, col = "blue", lty = 3)
lines(p_values_gamma, L3_ci_upper_gamma, lwd = 1.2, col = "blue", lty = 3)

# Legend
legend("topleft",
       c("Empirical L3", "True L3", "Bootstrap mean", "95% CI bounds", "Confidence area"),
       col = c("red", "purple", "green", "blue", rgb(0.3, 0.6, 1, 0.25)),
       lty = c(1, 1, 2, 3, NA), lwd = c(2.5, 2.5, 2, 1.2, NA),
       fill = c(NA, NA, NA, NA, rgb(0.3, 0.6, 1, 0.25)),
       bty = "o", cex = 0.95)

grid(col = "gray80", lty = 2)

dev.off()
cat("✓ Plot saved: figure/task4_part1_gamma.png\n\n")

# ============================================================================
# PART 2: NBA Salary Data by Marital Status
# ============================================================================

cat("=== PART 2: nbasal Wage Data by Marital Status ===\n\n")

# Load data
data(nbasal)
nbasal_clean <- subset(nbasal, complete.cases(wage, marr) & wage > 0)

cat("nbasal data loaded and cleaned\n")
cat("Total observations:", nrow(nbasal_clean), "\n")
cat("Married (marr=1):", sum(nbasal_clean$marr == 1), "\n")
cat("Unmarried (marr=0):", sum(nbasal_clean$marr == 0), "\n\n")

# Separate by marital status
married_wage <- nbasal_clean$wage[nbasal_clean$marr == 1]
unmarried_wage <- nbasal_clean$wage[nbasal_clean$marr == 0]

cat("Married wage - mean:", mean(married_wage), ", SD:", sd(married_wage), "\n")
cat("Unmarried wage - mean:", mean(unmarried_wage), ", SD:", sd(unmarried_wage), "\n\n")

# Evaluation points - same for both groups
p_values_wage <- seq(0.05, 0.95, length.out = 100)

# ---- Married group ----
cat("Bootstrap for married group (R=10000)...\n")
L3_original_married <- compute_l3(married_wage, p_values_wage)

boot_l3_married <- boot(
  data = married_wage,
  statistic = l3_bootstrap_stat,
  R = 10000,
  p_values = p_values_wage
)

L3_mean_married <- apply(boot_l3_married$t, 2, mean)
L3_ci_lower_married <- apply(boot_l3_married$t, 2, quantile, probs = 0.025)
L3_ci_upper_married <- apply(boot_l3_married$t, 2, quantile, probs = 0.975)

cat("✓ Married bootstrap complete\n")
cat("L3 range - original: [", min(L3_original_married), ", ", max(L3_original_married), "]\n")
cat("L3 range - CI: [", min(L3_ci_lower_married), ", ", max(L3_ci_upper_married), "]\n\n")

# ---- Unmarried group ----
cat("Bootstrap for unmarried group (R=10000)...\n")
L3_original_unmarried <- compute_l3(unmarried_wage, p_values_wage)

boot_l3_unmarried <- boot(
  data = unmarried_wage,
  statistic = l3_bootstrap_stat,
  R = 10000,
  p_values = p_values_wage
)

L3_mean_unmarried <- apply(boot_l3_unmarried$t, 2, mean)
L3_ci_lower_unmarried <- apply(boot_l3_unmarried$t, 2, quantile, probs = 0.025)
L3_ci_upper_unmarried <- apply(boot_l3_unmarried$t, 2, quantile, probs = 0.975)

cat("✓ Unmarried bootstrap complete\n")
cat("L3 range - original: [", min(L3_original_unmarried), ", ", max(L3_original_unmarried), "]\n")
cat("L3 range - CI: [", min(L3_ci_lower_unmarried), ", ", max(L3_ci_upper_unmarried), "]\n\n")

# ---- Overlap analysis ----
cat("=== Overlap Analysis ===\n\n")

# Check if confidence areas overlap
overlap_count <- sum(
  (L3_ci_upper_married >= L3_ci_lower_unmarried) & 
  (L3_ci_lower_married <= L3_ci_upper_unmarried)
)
overlap_pct <- 100 * overlap_count / length(p_values_wage)

cat("Number of p values where CIs overlap:", overlap_count, "out of", length(p_values_wage), "\n")
cat("Percentage of overlap: ", round(overlap_pct, 2), "%\n\n")

# Distance between curves
distance_original <- mean(abs(L3_original_married - L3_original_unmarried))
distance_mean <- mean(abs(L3_mean_married - L3_mean_unmarried))

cat("Mean absolute difference (original curves):", round(distance_original, 4), "\n")
cat("Mean absolute difference (bootstrap means):", round(distance_mean, 4), "\n\n")

# Visualization - Part 2 (combined plot)
png(file = "figure/task4_part2_wage.png", width = 10, height = 6, units = "in", res = 150)

y_min <- min(L3_ci_lower_married, L3_ci_lower_unmarried) - 0.05
y_max <- max(L3_ci_upper_married, L3_ci_upper_unmarried) + 0.05

# Empty plot first to avoid layering issues
plot(p_values_wage, L3_original_married,
     type = "n",
     main = "L3 Lorenz Curves by Marital Status\nwith 95% Bootstrap Confidence Areas",
     xlab = "p", ylab = "L3(p; Q)",
     ylim = c(y_min, y_max),
     cex.main = 1.2, cex.lab = 1.1)

# Confidence area - married (light red)
polygon(
  c(p_values_wage, rev(p_values_wage)),
  c(L3_ci_lower_married, rev(L3_ci_upper_married)),
  col = rgb(1, 0.5, 0.5, alpha = 0.25),
  border = NA
)

# Confidence area - unmarried (light blue)
polygon(
  c(p_values_wage, rev(p_values_wage)),
  c(L3_ci_lower_unmarried, rev(L3_ci_upper_unmarried)),
  col = rgb(0.3, 0.6, 1, alpha = 0.25),
  border = NA
)

# Original empirical L3 curves
lines(p_values_wage, L3_original_married, lwd = 2.5, col = "darkred", lty = 1)
lines(p_values_wage, L3_original_unmarried, lwd = 2.5, col = "darkblue", lty = 1)

# Bootstrap mean curves
lines(p_values_wage, L3_mean_married, lwd = 2, col = "darkred", lty = 2)
lines(p_values_wage, L3_mean_unmarried, lwd = 2, col = "darkblue", lty = 2)

# CI bounds
lines(p_values_wage, L3_ci_lower_married, lwd = 1.2, col = "red", lty = 3)
lines(p_values_wage, L3_ci_upper_married, lwd = 1.2, col = "red", lty = 3)
lines(p_values_wage, L3_ci_lower_unmarried, lwd = 1.2, col = "blue", lty = 3)
lines(p_values_wage, L3_ci_upper_unmarried, lwd = 1.2, col = "blue", lty = 3)

# Legend
legend("topleft",
       c("Married - Original", "Unmarried - Original",
         "Married - Bootstrap mean", "Unmarried - Bootstrap mean",
         "CI bounds", "Confidence areas"),
       col = c("darkred", "darkblue", "darkred", "darkblue", "gray50", "gray50"),
       lty = c(1, 1, 2, 2, 3, NA), lwd = c(2.5, 2.5, 2, 2, 1.2, NA),
       fill = c(NA, NA, NA, NA, NA, "gray50"),
       bty = "o", cex = 0.9)

grid(col = "gray80", lty = 2)

dev.off()
cat("✓ Plot saved: figure/task4_part2_wage.png\n\n")

# ============================================================================
# Save results
# ============================================================================

cat("=== Summary of Key Findings ===\n\n")

cat("PART 1 - Gamma(α=2, β=4):\n")
cat("  Original L3 range: [", round(min(L3_original_gamma), 4), ", ", round(max(L3_original_gamma), 4), "]\n")
cat("  Bootstrap mean range: [", round(min(L3_mean_gamma), 4), ", ", round(max(L3_mean_gamma), 4), "]\n")
cat("  Average bias: ", round(mean(L3_mean_gamma - L3_original_gamma), 4), "\n\n")

cat("PART 2 - NBA Salary by Marital Status:\n")
cat("  Married L3 range: [", round(min(L3_original_married), 4), ", ", round(max(L3_original_married), 4), "]\n")
cat("  Unmarried L3 range: [", round(min(L3_original_unmarried), 4), ", ", round(max(L3_original_unmarried), 4), "]\n")
cat("  Mean distance between curves: ", round(distance_original, 4), "\n")
cat("  CI overlap: ", overlap_pct, "% of p values\n\n")

cat("  Interpretation:\n")
if (overlap_pct >= 80) {
  cat("  → Confidence areas substantially overlap (descriptively suggesting similar curves).\n")
} else if (overlap_pct >= 50) {
  cat("  → Confidence areas partially overlap (descriptively suggesting some difference).\n")
} else {
  cat("  → Confidence areas show limited overlap (descriptively suggesting substantive difference).\n")
}
cat("  Note: Pointwise confidence areas are descriptive and do not form a formal global test.\n")

cat("\n")

# Save complete results
task4_results <- list(
  # Part 1 - Gamma
  gamma = list(
    p_values = p_values_gamma,
    L3_original = L3_original_gamma,
    L3_mean = L3_mean_gamma,
    L3_ci_lower = L3_ci_lower_gamma,
    L3_ci_upper = L3_ci_upper_gamma,
    boot_object = boot_l3_gamma,
    sample = gamma_sample,
    sample_size = n_gamma,
    alpha_true = alpha_true,
    beta_true = beta_true
  ),
  # Part 2 - Wage data
  wage = list(
    p_values = p_values_wage,
    married = list(
      n = length(married_wage),
      L3_original = L3_original_married,
      L3_mean = L3_mean_married,
      L3_ci_lower = L3_ci_lower_married,
      L3_ci_upper = L3_ci_upper_married,
      boot_object = boot_l3_married
    ),
    unmarried = list(
      n = length(unmarried_wage),
      L3_original = L3_original_unmarried,
      L3_mean = L3_mean_unmarried,
      L3_ci_lower = L3_ci_lower_unmarried,
      L3_ci_upper = L3_ci_upper_unmarried,
      boot_object = boot_l3_unmarried
    ),
    overlap_count = overlap_count,
    overlap_percentage = overlap_pct,
    mean_distance = distance_original
  )
)

saveRDS(task4_results, file = "results/task4_results.rds")
cat("✓ Results saved: results/task4_results.rds\n")
