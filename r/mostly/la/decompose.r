library(Matrix)

# create by 列向量
A <- t(matrix(c(3, -7, -2, 2, -3, 5, 1, 0, 6, -4, 0, -5, -9, 5, -5, 12), 4, 4))

# PA = LU
res <- Matrix::expand(Matrix::lu(A))
res$L

# Cholesky


# QR

# SVD
