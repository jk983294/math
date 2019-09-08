x <- c(9, 13, 21, 8, 36, 22, 12, 41, 31, 33, 19)
hist(x, xlab = "Weight")
hist(x, xlab = "Weight", col = "green", border = "red", xlim = c(0, 40), ylim = c(0, 
    5), breaks = 5)

# 轴须图
hist(x, xlab = "Weight")
rug(jitter(x))
lines(density(x), col = "blue", lwd = 2)

# 添加正态密度曲线和外框
h <- hist(x, xlab = "Weight")
xfit <- seq(min(x), max(x), length = 40)
yfit <- dnorm(xfit, mean = mean(x), sd = sd(x))
yfit <- yfit * diff(h$mids[1:2]) * length(x)
lines(xfit, yfit, col = "blue", lwd = 2)

# kernel density estimation
x <- c(5, 10, 15)
d <- density(x)
plot(d)
