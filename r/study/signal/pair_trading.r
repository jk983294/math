source("utils/utils.r")

run_regression <- function(dF) {
    return(coef(lm(y ~ x - 1, data = as.data.frame(dF))))
}
rolling_beta <- function(z, width) {
    rollapply(z, width = width, FUN = run_regression, by.column = FALSE, align = "right")
}

## TODO price rehabilitation
dat1 <- read_stock_data("../data/000001.csv")
dat2 <- read_stock_data("../data/600848.csv")
merged_close <- merge(dat1$close, dat2$close)  # two return matrix
filtered_close <- remove_missing_rows(merged_close)
x <- filtered_close$close
y <- filtered_close$close.1
start_date <- "2017-05-01"
end_date <- "2018-09-01"
start_date_out_sample <- "2018-09-01"
end_date_out_sample <- "2019-09-01"
range <- paste(start_date, "::", end_date, sep = "")
range_out_sample <- paste(start_date_out_sample, "::", end_date_out_sample, sep = "")

## seems static spread not promising, try dynamic beta Rolling window of trading days
window_length <- 10
x <- filtered_close$close[range]
y <- filtered_close$close.1[range]
dF <- cbind(x, y)
names(dF) <- c("x", "y")
betas <- rolling_beta(diff(dF), window_length)
data <- merge(betas, dF)
data$spread <- data$y - lag(betas, 1) * data$x
returns <- diff(dF)/dF
return_beta <- rolling_beta(returns, window_length)
data$spreadR <- diff(data$y)/data$y - return_beta * diff(data$x)/data$x
threshold <- sd(data$spread, na.rm = TRUE)
plot(data$spread, main = "AAPL vs. SPY In-Sample", cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)
abline(h = threshold, lty = 3)
abline(h = -threshold, lty = 2)

# out sample
x <- filtered_close$close[range_out_sample]
y <- filtered_close$close.1[range_out_sample]
dF <- cbind(x, y)
names(dF) <- c("x", "y")
beta_out_of_sample <- rolling_beta(diff(dF), 10)
# Buy and sell threshold
data_out <- merge(beta_out_of_sample, dF)
data_out$spread <- data_out$y - lag(beta_out_of_sample, 1) * data_out$x
# Plot the spread with in-sample bands
plot(data_out$spread, main = "AAPL vs. SPY out of sample", cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)
abline(h = threshold, lwd = 2)
abline(h = -threshold, lwd = 2)
# Generate sell and buy signals
buys <- ifelse(data_out$spread > threshold, 1, 0)
sells <- ifelse(data_out$spread < -threshold, -1, 0)
data_out$signal <- buys + sells
# plot buy sell signals
plot(data_out$spread, main = "AAPL vs. SPY out of sample", cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)
abline(h = threshold, lty = 2)
abline(h = -threshold, lty = 2)
point_type <- rep(NA, nrow(data_out))
buy_index <- which(data_out$signal == 1)
sell_index <- which(data_out$signal == -1)
point_type[buy_index] <- 21
point_type[sell_index] <- 24
points(data_out$spread, pch = point_type)
# trading opportunity count
num_of_buy_signals <- sum(buys, na.rm = TRUE)
num_of_sell_signals <- sum(abs(sells), na.rm = TRUE)
