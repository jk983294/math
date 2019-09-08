a <- matrix(c(1, 2, 3, 4), nrow = 2)  # [[1, 3], [2, 4]]
b <- matrix(c(1, -1, 0, 1), nrow = 2)  # [[1, 0], [-1, 1]]
x <- c(2, 5)

2 * a  # 2 * a = [[2, 6], [4, 8]]
t(a)  # transpose [[1, 2], [3, 4]]
a + b  # a + b = [[2, 3], [1, 5]]
a %*% x  # a * x = [[17], [24]]
a %*% b  # a * b = [[-2, 3], [-2, 4]]
rowMeans(a)  # [2 3]
rowSums(a)  # [4 6]
colMeans(a)  # [1.5 3.5]
colSums(a)  # [3 7]
crossprod(a, b)  # t(a) %*% b = [[-1, 2], [-1, 4]]
tcrossprod(a, b)  # a %*% t(b) = [[1, 2], [2, 2]]
solve(a)  # inverse of a = [[-2, 1.5], [1, -0.5]]
det(a)  # determinant -2
sweep(a, 1, c(1, 4), "+")  # # [[2, 4], [6, 8]] first row + 1, second row + 4

# special matrix/vector
matrix(0, nrow = 2, ncol = 2)  # 0-matrix [[0, 0], [0, 0]]
matrix(1, nrow = 2, ncol = 2)  # 1-matrix [[1, 1], [1, 1]]
diag(c(1, 2))  # diagonal matrix [[1, 0], [0, 2]]
diag(1, 2)  # 2*2 identity matrix [[1, 0], [0, 1]]
diag(2)  # 2*2 identity matrix [[1, 0], [0, 1]]
diag(a)  # diag vector [1 4]

# eigen
r <- eigen(a)
r$val  # eigen value [5.3723 -0.3723]
r$vec  # eigen vector

# svd
r <- svd(a)  # a = udv'
r$d  # vector containing the singular values
r$u  # matrix with columns contain the left singular vectors
r$v  # matrix with columns contain the right singular vectors

# cholesky decompose of a real symmetric positive-definite square matrix
r <- chol(matrix(c(1, 1, 1, 3), nrow = 2))  # upper triangular factor, such that r'r = a

# QR decomposition of a matrix
r <- qr(a)
r$qr  # an upper triangle that contains the decomposition and a lower triangle that contains information on the Q decomposition
r$rank  # the rank of A.
r$qraux  # a vector which contains additional information on Q
r$pivot  # contains info
