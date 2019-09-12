# period (N), initial value (x0), drift (mu), and variance
RW <- function(N, x0, mu, variance) {
    z <- cumsum(rnorm(n = N, mean = 0, sd = sqrt(variance)))
    t <- 1:N
    x <- x0 + t * mu + z
    return(x)
}
# mu is the drift

P1 <- RW(100, 10, 0, 4e-04)
P2 <- RW(100, 10, 0, 4e-04)
plot(P1, main = "Random Walk without Drift", xlab = "t", ylab = "Price", ylim = c(9.7, 10.3), typ = "l", col = "red")
par(new = T)  # to draw in the same plot
plot(P2, main = "Random Walk without Drift", xlab = "t", ylab = "Price", ylim = c(9.7, 10.3), typ = "l", col = "blue")
par(new = F)
