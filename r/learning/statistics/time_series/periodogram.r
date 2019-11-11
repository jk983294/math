library(TSA)

xk = sin(2 * pi * 0.2 * (0:99))
plot(xk, type = "h")
xkper <- periodogram(xk)
xkf <- fft(xk)
plot(xkf)
xkf[20:22]  # X[21] has big value, others small

# spectal leakage
xk1 = sin(2 * pi * 0.2 * (0:102))
xkper1 <- periodogram(xk1)  # period not complete
