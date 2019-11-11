set.seed(100)
X <- rnorm(1e+06, mean = 2.33, sd = 0.5)  # population
mu <- mean(X)
sd <- sd(X)
hist(X, breaks = 100)
abline(v = mu, lwd = 3, lty = 2)

# perform sampling from above population
mean_list <- list()
for (i in 1:10000) {
    # exercise such sample experiments
    mean_list[[i]] <- mean(sample(X, 10, replace = TRUE))  # sample from population
}
hist(unlist(mean_list), breaks = 500, xlab = "Mean of 10 samples from X", main = "Convergence of sample distribution", cex.main = 0.8)
abline(v = mu, lwd = 3, col = "white", lty = 2)
