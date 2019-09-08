library(compiler)

x <- runif(1e+06)
y <- runif(1e+06)
z <- vector(length = 1e+06)

system.time(z <- x + y) # vectorization is fastest

f <- function() for (i in 1:length(x)) {
    z[i] <<- x[i] + y[i]
}
system.time(f())

cf <- cmpfun(f)
system.time(cf())
