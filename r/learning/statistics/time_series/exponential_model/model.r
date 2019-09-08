# simple exponential smoothing
library(forecast)
fit <- HoltWinters(nhtemp, beta = FALSE, gamma = FALSE) # 单指数模型
plot(fit)
fit

forecast(fit, 1) # 一步向前预测

plot(forecast(fit, 1), xlab = "Year", ylab = expression(paste("Temperature (", degree * 
    F, ")", )), main = "New Haven Annual Mean Temperature")

accuracy(forecast(fit))


# exponential smoothing with level, slope, and seasonal components
fit <- HoltWinters(log(AirPassengers))
plot(fit)
fit

accuracy(forecast(fit))

pred <- forecast(fit, 5)
pred
plot(pred, main = "Forecast for Air Travel", ylab = "Log(AirPassengers)", xlab = "Time")
pred$mean <- exp(pred$mean) # 用原始尺度预测
pred$lower <- exp(pred$lower)
pred$upper <- exp(pred$upper)
p <- cbind(pred$mean, pred$lower, pred$upper)
dimnames(p)[[2]] <- c("mean", "Lo 80", "Lo 95", "Hi 80", "Hi 95")
p


# automatic exponential forecasting with ets()
library(forecast)
fit <- ets(JohnsonJohnson)
fit
plot(forecast(fit), main = "Johnson and Johnson Forecasts", ylab = "Quarterly Earnings (Dollars)", 
    xlab = "Time")
