############################### setup ###############################

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

  # return mean and confidence intervals
  rbind(mean_theta, ci_theta)
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
    theta_hat_bs[,i] <- bs[1,]

    # extract confidence intervals
    ci_alpha[,i] <- bs[2:3,1]
    ci_beta[,i] <- bs[2:3,2]
  }

  # calculate coverage probability
  coverage_alpha <- alpha > ci_alpha[1,] & alpha < ci_alpha[2,]
  coverage_beta <- beta > ci_beta[1,] & beta < ci_beta[2,]

  coverage <- paste0(round(100 * c(mean(coverage_alpha), mean(coverage_beta)), 2), "%")

  # absolute errors
  error_ml <- abs(theta_hat - theta)
  error_bs <- abs(theta_hat_bs - theta)

  # estimate closer to true theta
  closer <- ifelse(error_ml < error_bs, "ML",
                    ifelse(error_bs < error_ml, "BS", "SAME"))
  
  closer_frequency <- cbind(
    ML   = rowSums(closer == "ML"),
    BS   = rowSums(closer == "BS"),
    same = rowSums(closer == "same")
  )

  closer_str <- apply(closer_frequency, 1, paste, collapse = ":")

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
    closer_str,
    smaller_mean_err,
    coverage,
    row.names = c("alpha", "beta")
  ),
  bs_results = as.data.frame(theta_hat_bs)
)
  
}

######################### run the experiment ########################
results_summary <- list()
results_bs <- list()
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

    k <- k + 1
  }
}

all_results <- do.call(rbind, results_summary)

bs_results <- do.call(rbind, results_bs)

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