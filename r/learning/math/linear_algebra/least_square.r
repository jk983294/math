X <- cbind(c(1, 1, 1, 1, 1), c(1, 2, 3, 4, 5))
y <- c(3.7, 4.2, 4.9, 5.7, 6)

# y = a + bx beta = inv(X^T*X) * X^T * y
beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
beta_hat
