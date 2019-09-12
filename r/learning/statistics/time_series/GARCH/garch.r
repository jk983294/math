library(rugarch)

a <- runif(1000, min = 0, max = 100)  #random number

Spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(1, 1)), distribution.model = "std")
Garch <- ugarchfit(spec = Spec, data = a)
Garch

coefficient <- coef(Garch)
volatility <- sigma(Garch)
long.run.variance <- uncvariance(Garch)
