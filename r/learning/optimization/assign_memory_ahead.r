set.seed(1234)
k <- 1e+05
x <- rnorm(k)

y <- 0
system.time(for (i in 1:length(x)) y[i] <- x[i]^2) # for & no memory allocation

y <- numeric(k)
system.time(for (i in 1:k) y[i] <- x[i]^2) # for & memory allocation

y <- numeric(k)
system.time(y <- x^2) # vectorization & memory allocation
