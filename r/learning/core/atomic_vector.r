a <- 1  # actually a <- c(1)

x <- c(1, 2, 3, 4, 5, 6, 7, 8)
class(x) # "numeric"
print(x)
attr(x, "dim") <- c(2, 4) # change vector to matrix
print(x)
class(x) # "matrix"
attributes(x) # $dim [1] 2 4
attr(x, "dimnames") <- list(c("A1", "A2"), c("B1", "B2", "B3", "B4"))
print(x)
attr(x, "dim") <- NULL # change matrix back to vector
class(x) # "numeric"
print(x)
