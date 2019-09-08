# 存在季节性因素的时间序列数据(如月度数据、季度数据等)可以被分解为趋势因子、季节性因子和随机因子
plot(AirPassengers)
lAirPassengers <- log(AirPassengers)
plot(lAirPassengers, ylab = "log(AirPassengers)")
# stl函数只能处理相加模型,但这也不算一个多严重的限制,因为相乘模型总可以通过对数变换转换成相加模型
fit <- stl(lAirPassengers, s.window = "period")
plot(fit)
fit$time.series  # 每个观测值各分解项的值
exp(fit$time.series)  # 转化为原始尺度, anti-log


opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 1))
library(forecast)
monthplot(AirPassengers, xlab = "", ylab = "")  # 月度图, 每个月份组成的子序列
seasonplot(AirPassengers, year.labels = "TRUE", main = "")  # 季节图, 以年份为子序列
par(opar)
