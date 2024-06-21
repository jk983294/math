x <- c(2L, 0L, 3L)
n <- sum(x)
p <- c(1 / 3, 1 / 3, 1 / 3)
dmultinom(x, n, p) # PMF of P(X1 = 2, X2 = 0, X3 = 3)
rmultinom(10, n, p) # each column is a draw from Mult dist
