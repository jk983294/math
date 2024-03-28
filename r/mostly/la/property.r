library(MASS)

m <- matrix(c(1, 1, 1, 3), nrow = 2)
m1 <- matrix(c(1, 0, 0, 3), nrow = 2)

# det
det(m)
det(t(m)) # det(A) = det(transpose(A))
det(3 * m) # k^n * det(A) = det(k * A)
det(m %*% m1)
det(m) * det(m1) # det(AB) = det(A) * det(B)

# inverse
(inverse_m <- solve(m))

# Moore-Penrose inverse
# 1) A * inv(A) * A = A
# 2) inv(A) * A * inv(A) = inv(A)
# 3)


# eigen vector / value
