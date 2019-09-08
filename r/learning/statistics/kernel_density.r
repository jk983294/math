# default kernel is Gaussian
d <- density(mtcars$mpg)  # returns the density data
plot(d)  # plots the results

d <- density(mtcars$mpg)
plot(d, main = "Kernel Density of Miles Per Gallon")
polygon(d, col = "red", border = "blue")
rug(mtcars$mpg, col = "brown")

# comparing kernel density plots, draw them in same chart
set.seed(0)
x <- rnorm(100, 0, 1)
y <- rnorm(126, 0.3, 1.2)
z <- rnorm(93, -0.5, 0.7)
xlim <- range(x, y, z)
dx <- density(x, from = xlim[1], to = xlim[2], n = 200)
dy <- density(y, from = xlim[1], to = xlim[2], n = 200)
dz <- density(z, from = xlim[1], to = xlim[2], n = 200)

ylim <- range(dx$y, dy$y, dz$y)
plot(dx$x, dx$y, col = 1, lwd = 2, type = "l", xlim = xlim, ylim = ylim)
lines(dy$x, dy$y, col = 2, lwd = 2)
lines(dz$x, dz$y, col = 3, lwd = 2)

# using sm library
library(sm)
group.index <- rep(1:3, c(length(x), length(y), length(z)))
den <- sm.density.compare(c(x, y, z), group = group.index, model = "equal")
den2 <- sm.density.compare(c(x, y), group = rep(1:2, c(length(x), length(y))), model = "equal")
