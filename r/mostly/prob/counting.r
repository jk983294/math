A <- function(n, k) {
    factorial(n) / factorial(n - k)
}
C <- function(n, k) {
    choose(n, k)
}

factorial(9) # 9!
lfactorial(9) # log factorial for big number
choose(9, 2) # C(n, k)
lchoose(9, 2) # log C(n, k)
C(9, 2) # C(n, k)
A(9, 2) # A(n, k)

# x_1 + ... + x_r = n 的非负整数解的个数 C(n + r - 1, r - 1)
n <- 4
r <- 2
(ans <- choose(n + r - 1, r - 1))


combn(4, 2) # list all combination of C(n, k)
combinat::combn(4, 2) # list all combination of C(n, k)
combinat::permn(3) # list all combination of n!

