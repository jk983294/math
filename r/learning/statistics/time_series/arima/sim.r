library(forecast)
library(TSA)
library(tseries)

# sim ar2 model
ar.sim <- arima.sim(model = list(ar = c(0.9, -0.2)), n = 1000)
ts.plot(ar.sim)
ar.acf <- acf(ar.sim, type = "correlation", plot = T)
ar.pacf <- acf(ar.sim, type = "partial", plot = T)


# sim ma2 model
ma.sim <- arima.sim(model = list(ma = c(-0.7, 0.1)), n = 1000)
ts.plot(ma.sim)
ma.acf <- acf(ma.sim, type = "correlation", plot = T)
ma.pacf <- acf(ma.sim, type = "partial", plot = T)

# sim arma22 model
arma.sim <- arima.sim(model = list(ar = c(0.9, -0.2), ma = c(-0.7, 0.1)), n = 1000)
ts.plot(arma.sim)
arma.acf <- acf(arma.sim, type = "correlation", plot = T)
arma.pacf <- acf(arma.sim, type = "partial", plot = T)

# fit arma model
arma22mod <- arima(x = arma.sim, order = c(2, 0, 2))
tsdiag(arma22mod)
arma22mod$aic  # different fit parameters, choose min aic


# trend
data("co2", package = "TSA")
plot(co2)
tsdisplay(co2)
plot(diff(co2))  # remove trend component
attributes(co2)
structure(co2)  # check frequency of underling data
fit <- tslm(co2 ~ trend + season)
attributes(fit)
fit$coefficients
plot(fit$residuals)
periodogram(fit$residuals)  # check period frequency
plot(stl(co2, s.window = "per"))  # decompose trend season

# play with stock index
plot(EuStockMarkets[, "SMI"])

# high order trend fit
data("co2", package = "datasets")
zk <- co2
plot(zk)
tsdisplay(zk, plot.type = "partial")
tsdisplay(zk, plot.type = "spectrum")
tvec <- 1:NROW(zk)  # time vector
trfit <- lm(zk ~ tvec)  # poly fit first order
trfit2 <- lm(zk ~ I(tvec) + I(tvec^2))  # poly fit second order
trfit3 <- lm(zk ~ I(tvec) + I(tvec^2) + I(tvec^3))  # poly fit third order
summary(trfit)
plot(ts(trfit$residuals, start = tsp(zk)[1], frequency = 12), ylab = "residual")
periodogram(ts(trfit$residuals, start = tsp(zk)[1], frequency = 12))
plot(ts(trfit2$residuals, start = tsp(zk)[1], frequency = 12), ylab = "residual")
periodogram(ts(trfit2$residuals, start = tsp(zk)[1], frequency = 12))
plot(ts(trfit3$residuals, start = tsp(zk)[1], frequency = 12), ylab = "residual")
periodogram(ts(trfit3$residuals, start = tsp(zk)[1], frequency = 12))

plot(decompose(zk, type = "additive"))  # additive model, v[k] = trend[k] + season[k] + w[k], multiply mode, v[k] = trend[k] * season[k] * w[k]
plot(decompose(EuStockMarkets[, "SMI"], type = "additive"))  # random is not stationary, something wrong

# integration process
vk <- arima.sim(model = list(order = c(0, 1, 0)), n = 1000)
var(vk)
periodogram(vk)  # low frequency dominate
var(diff(vk))  # should << var(vk)
periodogram(diff(vk))
adf.test(vk)  # test unit root 1
adf.test(diff(vk))
