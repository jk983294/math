# transforming the time series and assessing stationarity
library(forecast)
library(tseries)
plot(Nile)
ndiffs(Nile)  # 1, exist order one trend
dNile <- diff(Nile)  # 差分一次, 下降的趋势被移除了
plot(dNile)
adf.test(dNile)  # 检验结果显示序列此时是平稳的
Acf(dNile) # 自相关图
Pacf(dNile) # 偏自相关图

# fit an ARIMA model, params(0, 1, 1) from Acf and Pacf
fit <- arima(Nile, order = c(0, 1, 1))
fit # AIC 值越小越好
accuracy(fit)


# evaluating the model fit
qqnorm(fit$residuals) # 数据满足正态分布,则数据中的点会落在图中的线上
qqline(fit$residuals)
Box.test(fit$residuals, type = "Ljung-Box") # 检验残差的自相关系数是否都为零


# forecasting with an ARIMA model
forecast(fit, 3)
plot(forecast(fit, 3), xlab = "Year", ylab = "Annual Flow")


# automated ARIMA forecasting
library(forecast)
fit <- auto.arima(sunspots)
fit
forecast(fit, 3)
accuracy(fit)
