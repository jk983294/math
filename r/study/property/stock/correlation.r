source("utils/utils.r")

dat1 <- read_stock_data("../data/000001.csv")
returns1 <- get_log_return(dat1$close)
returns2 <- get_log_return(dat1$open)
sv <- merge(returns1, returns2)  # two return matrix
sv <- remove_missing_rows(sv)
outliers <- which(sv[, 2] > 1)  # find the outliers
if (length(outliers) > 0) {
    sv <- sv[-outliers, ]  # if any outliers exist, remove them
}
cor(sv)  # open close correlation

## regression, use open to predict close
model <- lm(close ~ open, data = sv)
print(model)
print(summary(model))  # error, residual etc

par(mfrow = c(2, 2))
plot(model$residuals, main = "Residuals through time", xlab = "Days", ylab = "Residuals")  # scatter of residuals
hist(model$residuals, breaks = 100, main = "Distribution of residuals", xlab = "Residuals")
qqnorm(model$residuals)
qqline(model$residuals)
acf(model$residuals, main = "Autocorrelation")
par(mfrow = c(1, 1))

## check if one lead another, cross correlation
ccf(as.numeric(sv[, 1]), as.numeric(sv[, 2]), main = "Cross correlation between open and close", ylab = "Cross correlation", xlab = "Lag", 
    cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## volatility has memory effect
volatility <- sv[, 1]^2  # approximate of stddev of return since mean of return is 0
acf(volatility, main = "Actual returns squared", cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)
acf(sv[, 1]^3)  # L3 norm of return
acf(abs(sv[, 1]))  # L1 norm of return
