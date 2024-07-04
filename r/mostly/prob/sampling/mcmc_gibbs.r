# Define the parameters of the bivariate normal distribution
cov_ <- matrix(c(1, 0.5, 0.5, 1), nrow = 2) # covariance matrix

mu1 <- 0.
mu2 <- 0.
sigma1 <- cov_[1, 1]
sigma2 <- cov_[2, 2]
rho <- cov_[1, 2]
n_samples <- 10000
samples <- matrix(0, nrow = n_samples, ncol = 2)

# Initialize the first sample
x2 <- rnorm(1, mu2, sqrt(sigma2))

for (i in 1:n_samples) {
    # Sample x1 given x2
    mean_x1 <- mu1 + rho / sigma2 * (x2 - mu2)
    var_x1 <- sigma1 - rho^2 / sigma2
    x1 <- rnorm(1, mean_x1, sqrt(var_x1))

    # Sample x2 given x1
    mean_x2 <- mu2 + rho / sigma1 * (x1 - mu1)
    var_x2 <- sigma2 - rho^2 / sigma1
    x2 <- rnorm(1, mean_x2, sqrt(var_x2))

    # Store the samples
    samples[i, ] <- c(x1, x2)
}

# Burn-in period: discard the first 1000 samples
burn_in <- 1000
samples <- samples[(burn_in + 1):n_samples, ]

# Estimate the mean and covariance matrix from the samples
(estimated_mean <- colMeans(samples))
(estimated_covariance <- cov(samples))
