# theoretical acf/pacf
acfma1 <- ARMAacf(ma = 0.6, lag.max = 10)
plot((0:10), acfma1, type = "h")

pacfma1 <- ARMAacf(ma = 0.6, lag.max = 10, pacf = TRUE)
plot((0:9), pacfma1, type = "h")

acfar1 <- ARMAacf(ar = 0.6, lag.max = 10)
plot((0:10), acfar1, type = "h")

pacfar1 <- ARMAacf(ar = 0.6, lag.max = 10, pacf = TRUE)
plot((0:9), pacfar1, type = "h")

# estimate white noise acf/pacf is a impulse at 0
ek = rnorm(1000)
acf(ek)
pacf(ek)
