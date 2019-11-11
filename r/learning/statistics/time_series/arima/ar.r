# AR model
ar_coeff <- c(0.6, 0.3)
order <- length(ar_coeff)
ek = rnorm(1000)
# alpha = 0.6
vk <- 1:1000
for (i in 1:order) {
    vk[i] <- ek[i]
}
for (i in (order + 1):1000) {
    vk[i] <- ek[i]
    for (j in 1:order) {
        vk[i] <- vk[i] + vk[i - j] * ar_coeff[j]
    }
}
plot(vk, type = "l")  # check stationarity
points(vk, pch = "*")
summary(vk)
hist(vk)  # check gaussianity
acf(vk)
pacf(vk)

# fit to known order
armodel <- ar(vk, aic = FALSE, order.max = order)
armodel$resid[1:5]  # missing first p observation
acf(armodel$resid, na.action = na.pass)  # check residule is white noise
sigmaT <- armodel$asy.var.coef  # order*order matrix
sqrt(diag(sigmaT))  # parameter estimate error
armodel$ar  # parameter estimated

# find unknown order
m1 <- ar(vk, method = "mle")
m1$order
m2 <- arima(vk, order = c(m1$order, 0, 0))
sqrt(m2$sigma2)  # residual std err
p1 <- c(1, -m2$coef[1:2])  # characteristic equation
roots <- polyroot(p1)  # find solution
Mod(roots)  # absolute value of solution
