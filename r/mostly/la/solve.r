X <- rbind(
    c(1., -2, 1),
    c(0, 2, -8),
    c(-4, 5, 9)
)
y <- c(0., 8, 9)

beta <- solve(X, y)
inverse_m <- solve(X)

# backsolve (R, b) 求解上三角方程组 Rx = b
backsolve(matrix(c(1, 0, 2, 2.), 2, 2), c(1, 2))

# forwardsolve(L, b) 求解下三角方程组 Lx = b
forwardsolve(matrix(c(1, 2, 0, 2.), 2, 2), c(1, 3))

library(matlib)
# 有唯一的解, the solution is unique if r(A|b)=r(A)=n
matlib::showEqn(X, y)
matlib::echelon(X, y)
matlib::Solve(X, y, fractions = TRUE)


# 有无数的解, the solution is underdetermined if r(A|b)=r(A)<n
A <- matrix(c(1, -1, -2, 2), 2, 2)
b <- c(-1, 1)
matlib::showEqn(A, b)
matlib::echelon(A, b)
matlib::Solve(A, b, fractions = TRUE)
matlib::plotEqn(A, b, xlim = c(-2, 4))

# 无解, the equations are inconsistent (no solutions) if r(A|b)>r(A)
A <- matrix(c(1, 2, 3, -1, 2, 1), 3, 2)
b <- c(2, 1, 6)
matlib::showEqn(A, b)
matlib::echelon(A, b)
matlib::Solve(A, b, fractions = TRUE)
matlib::plotEqn(A, b, xlim = c(-2, 4))
