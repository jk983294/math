# used to estimate EX = ∑p(x)*f(x)
# μ = (∑w_i*x_i) / (∑w_i) where w_i = p(x_i) / q(x_i)
# q(x) is The distribution which we can easily sample
# we choose q(x) to be high where |p(x)*f(x)| is high
# this helps to reduce the variance of estimated mean, which means converge fast

# estimate the mean of a standard normal distribution using importance sampling with a uniform proposal distribution.
target_distribution <- function(x) { # target distribution (standard normal)
    return(dnorm(x, mean = 0, sd = 1))
}

proposal_distribution <- function(x) {
    return(dunif(x, min = -5, max = 5))
}

# Generate samples from the proposal distribution
n_samples <- 10000
samples <- runif(n_samples, min = -5, max = 5)

# Compute importance weights
weights <- target_distribution(samples) / proposal_distribution(samples)

# Normalize weights
normalized_weights <- weights / sum(weights)

# Estimate the mean of the target distribution
estimated_mean <- sum(normalized_weights * samples)

# Print the estimated mean
print(estimated_mean)
