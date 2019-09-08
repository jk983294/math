# Matrix Computations
matrix1 <- matrix(c(3, 9, -1, 4, 2, 6), nrow = 2)
matrix2 <- matrix(c(5, 2, 0, 9, 3, 4), nrow = 2)

# transpose
t(matrix1)

# element-wise computation
matrix1 + matrix2
matrix1 - matrix2
matrix1 * matrix2
matrix1/matrix2

# matrix multiply
M = matrix(c(2, 6, 5, 1, 10, 4), nrow = 2, ncol = 3, byrow = TRUE)
M %*% t(M)  # multiply with its transpose
2 * M  # scalar multiply
