y <- 1:3
x <- c(2, 8, 9)
t(y)  # transpose [1 2 3]^T
2 * y  # [2 4 6]
x + y  # [3 10 12]
sum(x * y)  # 45, list, inner product
crossprod(x, y)  # 45, 1*1 matrix, it actually is inner product
sqrt(sum(x * x))  # 12.21, length, norm
rep(0, 3)  # 0-vector [0 0 0]
rep(1, 3)  # 1-vector [1 1 1]
