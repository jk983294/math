source("utils/utils.r")

dat1 <- read_stock_data("../data/000001.csv")
dat2 <- read_stock_data("../data/600848.csv")
returns1 <- get_log_return(dat1$close)
returns2 <- get_log_return(dat2$close)
price_changes <- merge(returns1, returns2)  # two return matrix
price_changes <- remove_missing_rows(price_changes)

# prices <- merge(get_log_return(dat1$close), get_log_return(dat2$close)) price_changes <- apply(prices, 2, diff)
price_changes_1 <- as.numeric(price_changes[, 1])
price_changes_2 <- as.numeric(price_changes[, 2])
plot(price_changes_1, price_changes_2, xlab = "A price changes", ylab = "B price changes", main = "A vs. B", cex.main = 0.8, cex.lab = 0.8, 
    cex.axis = 0.8)
abline(lm(price_changes_1 ~ price_changes_2))
abline(lm(price_changes_2 ~ price_changes_1), lty = 2)
grid()

# Total least squares regression
r <- prcomp(~price_changes_1 + price_changes_2)
slope <- r$rotation[2, 1]/r$rotation[1, 1]
intercept <- r$center[2] - slope * r$center[1]
abline(a = intercept, b = slope, lty = 3)  # Show the first principal component on the plot
