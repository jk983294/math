set.seed(42)

# Number of data points
n <- 1000

# True parameters for the GMM
mu1 <- 0.
sigma1 <- 2.
mu2 <- 5.
sigma2 <- 2.

# Generate data
z <- sample(1:2, n, replace = TRUE, prob = c(0.5, 0.5))
data <- rep(0., n)
data[z == 1] <- rnorm(sum(z == 1), mean = mu1, sd = sigma1)
data[z == 2] <- rnorm(sum(z == 2), mean = mu2, sd = sigma2)

# Plot the data
plot(data, col = z, xlab = "X1", ylab = "X2", main = "Synthetic GMM Data")


# Initialize parameters
mu_est1 <- rnorm(1)
mu_est2 <- rnorm(1)
sd_est1 <- 1.
sd_est2 <- 1.
P_case1 <- 0.5
P_case2 <- 0.5

# Number of iterations
max_iter <- 100
tol <- 1e-6
log_likelihood <- rep(0, max_iter)

for (iter in 1:max_iter) {
    # E-step: Compute responsibilities
    P_x_given_case1 <- dnorm(data, mean = mu_est1, sd = sd_est1)
    P_x_given_case2 <- dnorm(data, mean = mu_est2, sd = sd_est2)
    total_P_x <- P_x_given_case1 * P_case1 + P_x_given_case2 * P_case2
    P_case1_given_x <- (P_x_given_case1 * P_case1) / total_P_x
    P_case2_given_x <- (P_x_given_case2 * P_case2) / total_P_x

    # M-step: Update parameters
    new_mu_est1 <- sum(P_case1_given_x * data) / sum(P_case1_given_x)
    new_sd_est1 <- sqrt(sum(P_case1_given_x * (data - new_mu_est1)^2) / sum(P_case1_given_x))
    new_mu_est2 <- sum(P_case2_given_x * data) / sum(P_case2_given_x)
    new_sd_est2 <- sqrt(sum(P_case2_given_x * (data - new_mu_est2)^2) / sum(P_case2_given_x))
    new_P_case1 <- mean(P_case1_given_x)
    new_P_case2 <- mean(P_case2_given_x)

    mu_est1 <- new_mu_est1
    sd_est1 <- new_sd_est1
    mu_est2 <- new_mu_est2
    sd_est2 <- new_sd_est2
    P_case1 <- new_P_case1
    P_case2 <- new_P_case2
    likelihood <- sum(log(P_x_given_case1) + log(P_x_given_case2))

    print(c(new_mu_est1, new_sd_est1, new_mu_est2, new_sd_est2, new_P_case1, new_P_case2, likelihood))

    # Compute log-likelihood
    log_likelihood[iter] <- likelihood

    # Check for convergence
    if (iter > 1 && abs(log_likelihood[iter] - log_likelihood[iter - 1]) < tol) {
        cat("Converged after", iter, "iterations\n")
        break
    }
}