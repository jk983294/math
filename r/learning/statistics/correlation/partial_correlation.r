library(ppcor)

w <- rnorm(200)
v <- rnorm(200)
z <- rnorm(200)
x <- 2 * z + 3 * w
y <- z + v

# correlation
cor(cbind(x, y))
pcor(cbind(x, y, z))$estimate  # under partial correlation, x & y correlation drops fast

# covariance
cov(cbind(x, y))
cov2cor(cov(cbind(x, y)))
