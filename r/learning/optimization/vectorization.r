set.seed(1234)
mymatrix <- matrix(rnorm(1e+07), ncol = 10)
system.time(colSums(mymatrix))  # using vectorization
