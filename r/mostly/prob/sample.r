n <- 10L
k <- 5L
sample(n, k)
sample(n, k, replace = TRUE)
sample(-10:10, k, replace = TRUE)
sample(4, k, replace = TRUE, prob = c(0.1, 0.2, 0.3, 0.4)) # weighted sample
times <- 10^4
replicate(times, sum(sample(-10:10, k, replace = TRUE) > 0)) # simulation with times

# n个人至少2人生日相同的概率
# birthday simulation
b <- sample(1:365, 23, replace = TRUE)
tabulate(b)
r <- replicate(times, max(tabulate(sample(1:365, 23, replace = TRUE))))
(prob <- sum(r >= 2) / times)
pbirthday(23) # prob at least one match if n people
qbirthday(0.5) # number of people needed in order to have prob p

# match simulation
n <- 100L
r <- replicate(times, sum(sample(n) == (1:n)))
mean(r) # = 1.
