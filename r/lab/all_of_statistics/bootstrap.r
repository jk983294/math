library(moments)

n <- 100
sample_y <- rnorm(n)
sample_x <- exp(sample_y)
mean_hat <- mean(sample_x)
sd_hat <- sd(sample_x)
skew_hat <- sum((sample_x - mean_hat)^3)/n/(sd_hat^3) - 3

B <- 1000
Tboot <- 1:B

for (i in 1:B) {
    xx1 <- sample(sample_x, n, replace = TRUE)
    mean_hat_xx1 <- mean(xx1)
    sd_hat_xx1 <- sd(xx1)
    Tboot[i] <- sum((xx1 - mean_hat_xx1)^3)/n/(sd_hat_xx1^3) - 3
}

se <- sqrt(var(Tboot))
Normal.intvl <- c(skew_hat - 2 * se, skew_hat + 2 * se)
percentile.intvl <- c(quantile(Tboot, 0.025), quantile(Tboot, 0.975))
pivotal.intvl <- c(2 * skew_hat - quantile(Tboot, 0.975), 2 * skew_hat - quantile(Tboot, 
    0.025))
