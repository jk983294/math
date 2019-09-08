# a*x = b
a <- matrix(c(1, 1, -1, 1), nrow = 2)
b <- c(2, 4)
solve(a, b)  # [3 1] list
solve(a) %*% b  # [[3], [1]] matrix
