result <- array(c(5, 9, 3, 10, 11, 12, 13, 14, 15), dim = c(3, 3, 2))

# create matrices from these arrays.
matrix1 <- result[, , 1]
matrix2 <- result[, , 2]
matrix1 + matrix2
