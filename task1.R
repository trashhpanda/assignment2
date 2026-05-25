############################### setup ###############################

# Create output directories if they don't exist
dir.create("results", showWarnings = FALSE)
dir.create("figure", showWarnings = FALSE)

# set random seed
set.seed(42)

# number of MC iterations
MC <- 1000

# number of bootstrap iterations
BS <- 10000

# sample sizes
N <- c(20, 50, 200)

# true parameter values for the Gamma distribution 
theta_true <- matrix(
  c(0.5, 1,
    10, 1,
    10, 3),
  ncol = 2,
  byrow = TRUE
)


#################### maximum likelihood estimator ###################
# probably a bit overengineered
mle_theta <- function(x) {
  # calculate c = log(mean(xi)) - mean(log(xi)) part for all bootstrap samples
  xbar = colMeans(x)
  c <- log(xbar) - colMeans(log(x))

  # placeholder for alpha estimates
  alpha_hat <- numeric(length(c))

  # now find alpha for each bootstrap sample
  for (i in 1:length(c)){ 
    # score function
    f <- function(alpha) {
      log(alpha) - digamma(alpha) - c[i]
    }
  
    # find root of the score function
    alpha_hat[i] <- uniroot(
      f,
      interval = c(1e-6, 1000)
    )$root
  }
  
  # calculate beta estimates
  beta_hat <- xbar/alpha_hat

  matrix(c(alpha_hat, beta_hat), nrow = 2, byrow = TRUE)
}

############################# bootstrap #############################
bootstrap <- function(x, B) {
  # sample size
  n <- length(x)

  # draw n indices for B bootstrap samples
  idx <- matrix(
    sample.int(n, n*B, replace=TRUE),
    nrow = n
  )

  # B bootstrap samples of n observations
  bs_samples <- matrix(x[idx], nrow = n)

  # calculate mle for the bootstrap samples
  theta_hat <- mle_theta(bs_samples)

  # construct 95% confidence interval from all bootstrap estimates
  ci_theta <- apply(theta_hat, 1, quantile,
              probs = c(0.025, 0.975),
              type = 7)
  
  # calculate mean of the bootstrap estimates
  mean_theta <- rowMeans(theta_hat)

  # return full bootstrap estimates, mean, and confidence intervals
  list(
    theta_boot = theta_hat,
    mean_theta = mean_theta,
    ci_theta = ci_theta
  )
}

####################### monte carlo simulation ######################
mc_simulation <- function(n, theta, M, B) {
  # get shape and scale
  alpha <- theta[1]
  beta <- theta[2]

  # mle results
  theta_hat <- matrix(NA, nrow=2, ncol=M)

  # mean bootstrap mle estimates
  theta_hat_bs <- matrix(NA, nrow=2, ncol=M)

  # bootstrap confidence intervals
  ci_alpha <- matrix(NA, nrow=2, ncol=M)
  ci_beta <- matrix(NA, nrow=2, ncol=M)
  
  # Store FULL bootstrap estimates for visualization (pick one representative sample)
  # We'll store bootstrap estimates from the first MC iteration
  full_bootstrap_alpha <- NULL
  full_bootstrap_beta <- NULL

  # draw M samples of size n
  X <- matrix(
    rgamma(n*M, shape=alpha, scale=beta),
    nrow = n,
    ncol = M
  )
  
  # calculate mle for the MC samples
  theta_hat <- mle_theta(X)

  # perform bootstrap
  for (i in 1:M) {
    bs <- bootstrap(X[,i], B)
    
    # extract mean of the bootstrap estimates for i-th MC sample
    theta_hat_bs[,i] <- bs$mean_theta

    # extract confidence intervals
    ci_alpha[,i] <- bs$ci_theta[,1]
    ci_beta[,i] <- bs$ci_theta[,2]
    
    # Store full bootstrap estimates from first iteration for visualization
    if (i == 1) {
      full_bootstrap_alpha <- bs$theta_boot[1, ]
      full_bootstrap_beta <- bs$theta_boot[2, ]
    }
  }

  # calculate coverage probability
  coverage_alpha <- alpha >= ci_alpha[1,] & alpha <= ci_alpha[2,]
  coverage_beta <- beta >= ci_beta[1,] & beta <= ci_beta[2,]

  coverage_alpha_prob <- mean(coverage_alpha)
  coverage_beta_prob <- mean(coverage_beta)

  # absolute errors
  error_ml <- abs(theta_hat - theta)
  error_bs <- abs(theta_hat_bs - theta)

  # estimate closer to true theta
  closer <- ifelse(error_ml < error_bs, "ML",
                    ifelse(error_bs < error_ml, "BS", "SAME"))
  
  closer_frequency <- cbind(
    ML   = rowSums(closer == "ML"),
    BS   = rowSums(closer == "BS"),
    SAME = rowSums(closer == "SAME")
  )

  # mean absolute error
  mean_err_ml <- rowMeans(error_ml)
  mean_err_bs <- rowMeans(error_bs)

  smaller_mean_err <- ifelse(mean_err_ml < mean_err_bs, "ML",
                      ifelse(mean_err_bs < mean_err_ml, "BS", "SAME"))
  
  # mean theta
  mean_theta_ml <- rowMeans(theta_hat)
  mean_theta_bs <- rowMeans(theta_hat_bs)

 
list(
  summary = data.frame(
    theta,
    mean_theta_ml,
    mean_theta_bs,
    mean_err_ml,
    mean_err_bs,
    closer_ML = closer_frequency[,"ML"],
    closer_BS = closer_frequency[,"BS"],
    closer_SAME = closer_frequency[,"SAME"],
    smaller_mean_err,
    coverage_alpha = coverage_alpha_prob,
    coverage_beta = coverage_beta_prob,
    row.names = c("alpha", "beta")
  ),
  bs_results = data.frame(
    alpha = theta_hat_bs[1, ],
    beta = theta_hat_bs[2, ]
  ),
  full_bootstrap_alpha = full_bootstrap_alpha,
  full_bootstrap_beta = full_bootstrap_beta
)
  
}

######################### run the experiment ########################
results_summary <- list()
results_bs <- list()
full_bootstrap_data <- list()
k <- 1

for (n in N) {
  print(paste("sample size:", n))

  for (i in 1:nrow(theta_true)){
    print(theta_true[i,])

    mc_results <- mc_simulation(n, theta_true[i,], MC, BS)
    
    results_summary[[k]] <- cbind(
      n = n,
      theta_id = i,
      param = c("alpha", "beta"),
      mc_results$summary
    )

    results_bs[[k]] <- cbind(
      n = n,
      theta_id = i,
      mc_results$bs_results
    )
    
    # Store full bootstrap data for visualization (from first MC iteration)
    full_bootstrap_data[[k]] <- data.frame(
      n = n,
      theta_id = i,
      alpha_bootstrap = mc_results$full_bootstrap_alpha,
      beta_bootstrap = mc_results$full_bootstrap_beta
    )

    k <- k + 1
  }
}

all_results <- do.call(rbind, results_summary)

bs_results <- do.call(rbind, results_bs)

full_bootstrap_df <- do.call(rbind, full_bootstrap_data)

rownames(all_results) <- NULL

all_results$n <- factor(all_results$n,
                        levels = c(20, 50, 200),
                        labels = c("n = 20", "n = 50", "n = 200"))

all_results$theta_id <- factor(all_results$theta_id,
                               levels = 1:3,
                               labels = c("theta1", "theta2", "theta3"))

bs_results$n <- factor(bs_results$n,
                        levels = c(20, 50, 200),
                        labels = c("n = 20", "n = 50", "n = 200"))

bs_results$theta_id <- factor(bs_results$theta_id,
                               levels = 1:3,
                               labels = c("theta1", "theta2", "theta3"))

print(all_results)

saveRDS(all_results, file = "results/task1_results.rds")

saveRDS(bs_results, file = "results/task1_bootstrap.rds")

saveRDS(full_bootstrap_df, file = "results/task1_full_bootstrap.rds")

######################### Visualizations: Histograms and QQ-plots #########################

# Load full bootstrap data
full_bootstrap_df <- readRDS("results/task1_full_bootstrap.rds")

# Get unique scenarios
scenarios <- unique(full_bootstrap_df[, c("n", "theta_id")])

# Create histograms and qq-plots for each combination
for (row_idx in 1:nrow(scenarios)) {
  n <- scenarios$n[row_idx]
  theta_id <- scenarios$theta_id[row_idx]
  
  boot_subset <- full_bootstrap_df[
    full_bootstrap_df$n == n & full_bootstrap_df$theta_id == theta_id,
  ]
  
  alpha_boot <- boot_subset$alpha_bootstrap
  beta_boot <- boot_subset$beta_bootstrap
  
  true_alpha <- theta_true[theta_id, 1]
  true_beta <- theta_true[theta_id, 2]
  
  # Create filename
  fig_name_alpha <- sprintf("figure/task1_hist_qq_alpha_n%d_theta%d.png", n, theta_id)
  fig_name_beta <- sprintf("figure/task1_hist_qq_beta_n%d_theta%d.png", n, theta_id)
  
  # Alpha: Histogram and QQ-plot
  png(filename = fig_name_alpha, width = 1000, height = 500)
  par(mfrow = c(1, 2))
  
  hist(alpha_boot,
       main = sprintf("Bootstrap Alpha Estimates\nn = %d, theta = (%.1f, %.1f)", 
                      n, true_alpha, true_beta),
       xlab = "Alpha",
       ylab = "Frequency",
       col = "lightblue",
       border = "darkblue",
       breaks = 30)
  
  abline(v = true_alpha, col = "red", lwd = 2, lty = 2)
  
  legend("topright", 
         legend = sprintf("True alpha = %.1f", true_alpha),
         col = "red", lty = 2, lwd = 2)
  
  qqnorm(alpha_boot,
         main = sprintf("QQ-plot Alpha\nn = %d, theta = (%.1f, %.1f)", 
                        n, true_alpha, true_beta))
  
  qqline(alpha_boot, col = "red", lwd = 2)
  
  par(mfrow = c(1, 1))
  dev.off()
  
  # Beta: Histogram and QQ-plot
  png(filename = fig_name_beta, width = 1000, height = 500)
  par(mfrow = c(1, 2))
  
  hist(beta_boot,
       main = sprintf("Bootstrap Beta Estimates\nn = %d, theta = (%.1f, %.1f)", 
                      n, true_alpha, true_beta),
       xlab = "Beta",
       ylab = "Frequency",
       col = "lightgreen",
       border = "darkgreen",
       breaks = 30)
  
  abline(v = true_beta, col = "red", lwd = 2, lty = 2)
  
  legend("topright",
         legend = sprintf("True beta = %.1f", true_beta),
         col = "red", lty = 2, lwd = 2)
  
  qqnorm(beta_boot,
         main = sprintf("QQ-plot Beta\nn = %d, theta = (%.1f, %.1f)", 
                        n, true_alpha, true_beta))
  
  qqline(beta_boot, col = "red", lwd = 2)
  
  par(mfrow = c(1, 1))
  dev.off()
  
  cat(sprintf("Saved: %s and %s\n", fig_name_alpha, fig_name_beta))
}

cat("\n=== Visualization complete ===\n")

######################### Normality tests for bootstrap distributions #########################

# Shapiro-Wilk test in R works only for sample sizes between 3 and 5000.
# Since each bootstrap distribution has B = 10,000 values, we test a random subset of 5000 values.
# The test is used as an additional numerical diagnostic together with histograms and QQ-plots.

set.seed(123)

normality_results <- list()
k_norm <- 1

scenarios <- unique(full_bootstrap_df[, c("n", "theta_id")])

for (row_idx in 1:nrow(scenarios)) {
  n_current <- scenarios$n[row_idx]
  theta_id_current <- scenarios$theta_id[row_idx]
  
  boot_subset <- full_bootstrap_df[
    full_bootstrap_df$n == n_current & full_bootstrap_df$theta_id == theta_id_current,
  ]
  
  alpha_boot <- boot_subset$alpha_bootstrap
  beta_boot <- boot_subset$beta_bootstrap
  
  # Take at most 5000 observations for Shapiro-Wilk test
  alpha_test_sample <- sample(alpha_boot, size = min(5000, length(alpha_boot)))
  beta_test_sample <- sample(beta_boot, size = min(5000, length(beta_boot)))
  
  shapiro_alpha <- shapiro.test(alpha_test_sample)
  shapiro_beta <- shapiro.test(beta_test_sample)
  
  normality_results[[k_norm]] <- data.frame(
    n = n_current,
    theta_id = theta_id_current,
    param = "alpha",
    shapiro_W = unname(shapiro_alpha$statistic),
    shapiro_p_value = shapiro_alpha$p.value,
    normality_conclusion = ifelse(
      shapiro_alpha$p.value > 0.05,
      "No evidence against normality",
      "Evidence against normality"
    )
  )
  
  k_norm <- k_norm + 1
  
  normality_results[[k_norm]] <- data.frame(
    n = n_current,
    theta_id = theta_id_current,
    param = "beta",
    shapiro_W = unname(shapiro_beta$statistic),
    shapiro_p_value = shapiro_beta$p.value,
    normality_conclusion = ifelse(
      shapiro_beta$p.value > 0.05,
      "No evidence against normality",
      "Evidence against normality"
    )
  )
  
  k_norm <- k_norm + 1
}

normality_results <- do.call(rbind, normality_results)

print(normality_results)

saveRDS(normality_results, file = "results/task1_normality_tests.rds")

cat("\n=== Normality test summary ===\n")
cat("Shapiro-Wilk tests were performed on random subsets of up to 5000 bootstrap estimates.\n")
cat("Interpretation rule: p-value > 0.05 means no evidence against normality; p-value <= 0.05 suggests evidence against normality.\n")
cat("The test results should be interpreted together with histograms and QQ-plots, because for large samples even small deviations from normality may be statistically significant.\n")