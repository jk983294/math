Q <- matrix(c(
    1 / 3, 1 / 3, 1 / 3, 0,
    0, 0, 1 / 2, 1 / 2,
    0, 1, 0, 0,
    1 / 2, 0, 0, 1 / 2
), nrow = 4, ncol = 4, byrow = TRUE)
Q2 <- Q %*% Q
Q3 <- Q2 %*% Q
Q4 <- Q3 %*% Q
Q5 <- Q4 %*% Q

# choose eigen vector whose eigen value = 1, then normalize it to stationary distribution
evs <- eigen(t(Q))

# simulation to get stationary distribution
M <- nrow(Q) # how many state
nsim <- 1e4
x <- rep(0, nsim) # pre-allocate space for results of simulation
x[1] <- sample(1:M, 1) # choose random initial state
for (i in 2:nsim) {
    x[i] <- sample(M, 1, prob = Q[x[i - 1], ]) # choose x[i]'s state based on Q[x[i - 1], ] probability
}
x <- x[-(1:(nsim / 2))] # remove first half burn-in simulation which is not stable
stationary_distribution <- table(x) / length(x)
