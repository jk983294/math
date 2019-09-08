set.seed(1234)
mymatrix <- matrix(rnorm(1e+07), ncol = 10)
accum <- function(x) {
    sums <- numeric(ncol(x))
    for (i in 1:ncol(x)) {
        for (j in 1:nrow(x)) {
            sums[i] <- sums[i] + x[j, i]
        }
    }
}

system.time(accum(mymatrix))  # using loops
system.time(colSums(mymatrix))  # using vectorization
