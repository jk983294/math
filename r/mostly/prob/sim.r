# simulating order statistics
order_stats <- replicate(10^5, sort(rnorm(10)))
x9 <- order_stats[9, ]
mean(x9)
sd(x9)
hist(x9)

# time until HH vs HT
paste(sample(c("H", "T"), 100, replace = TRUE), collapse = "") # one sequence of coin toss
r <- replicate(10^3, paste(sample(c("H", "T"), 100, replace = TRUE), collapse = ""))
t <- stringr::str_locate(r, "HH")
mean(t[, 2]) # average waiting time for HH
mean(stringr::str_locate(r, "HT")[, 2]) # average waiting time for HT

# Monte Carlo estimate of π
nsim <- 100000L
x <- runif(nsim, -1., 1.)
y <- runif(nsim, -1., 1.)
pi <- 4 * sum(x^2 + y^2 < 1) / nsim
