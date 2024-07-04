target_distribution <- function(x) { # target distribution (standard normal)
    return(dnorm(x, mean = 0, sd = 1))
}

# let's use a uniform distribution as the proposal distribution.
# We need to choose the bounds and the constant M such that p(x)≤Mq(x)
proposal_distribution <- function(x) {
    return(dunif(x, min = -5, max = 5))
}

# Constant M
M <- 1. / sqrt(2 * pi) * 10 # * 10 to offset proposal_distribution's density 0.1

sample_proposal <- function() { # Sample from the proposal distribution
    return(runif(1, min = -5, max = 5))
}

# Rejection sampling
rejection_sampling <- function(n_samples) {
    samples <- c()
    while (length(samples) < n_samples) {
        x <- sample_proposal()
        u <- runif(1, min = 0, max = 1)
        if (u < target_distribution(x) / (M * proposal_distribution(x))) {
            samples <- c(samples, x)
        }
    }
    return(samples)
}

# Generate samples
n_samples <- 10000
samples <- rejection_sampling(n_samples)

# Plot the histogram of the samples
hist(samples, breaks = 50, freq = FALSE, main = "Rejection Sampling", xlab = "x", col = "lightblue")
curve(dnorm(x, mean = 0, sd = 1), from = -5, to = 5, add = TRUE, col = "red", lwd = 2)
