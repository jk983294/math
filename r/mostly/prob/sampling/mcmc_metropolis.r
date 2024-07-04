stable_distribution <- c(0.4, 0.1, 0.2, 0.3) # define the desired distribution function
n_state <- length(stable_distribution)
n_samples <- 10000L
burn_in <- 1000L
thin <- 10L
chain <- numeric(n_samples * thin + burn_in)
chain[1] <- sample(1:n_state, 1) # start with a random state
for (i in 2:(n_samples * thin + burn_in)) {
    current <- chain[i - 1]
    proposal <- sample(1:n_state, 1)
    s_i <- stable_distribution[current]
    s_j <- stable_distribution[proposal]
    # assume p(j,i) = p(i,j) = 1/4, which cancel out each other
    acceptance_ratio <- min(s_j / s_i, 1.)
    if (runif(1) < acceptance_ratio) {
        chain[i] <- proposal
    } else {
        chain[i] <- current
    }
}

# Remove burn-in and thin the chain
samples <- chain[(burn_in + 1):length(chain)]
samples <- samples[seq(1, length(samples), by = thin)]
prob_dist <- table(samples) / length(samples)
barplot(prob_dist, main = "Bar Plot of Value Counts", xlab = "Values", ylab = "Frequency", col = "blue")
