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

calculate_spread <- function(x, y, beta) {
    return(y - beta * x)
}

calculate_beta_and_level <- function(x, y, start_date, end_date) {
    # calculate the beta and level given start and end dates
    require(xts)
    time_range <- paste(start_date, "::", end_date, sep = "")
    x <- x[time_range]
    y <- y[time_range]
    dx <- diff(x[time_range])
    dy <- diff(y[time_range])
    r <- prcomp(~dx + dy)
    beta <- r$rotation[2, 1]/r$rotation[1, 1]
    spread <- calculate_spread(x, y, beta)
    names(spread) <- "spread"
    level <- mean(spread, na.rm = TRUE)
    outL <- list()
    outL$spread <- spread
    outL$beta <- beta
    outL$level <- level
    return(outL)
}

calculate_buy_sell_signals <- function(spread, beta, level, lower_threshold, upper_threshold) {
    # calculate buy and sell signals with upper and lower threshold
    buy_signals <- ifelse(spread <= level - lower_threshold, 1, 0)
    sell_signals <- ifelse(spread >= level + upper_threshold, 1, 0)
    # bind these vectors into a matrix
    output <- cbind(spread, buy_signals, sell_signals)
    colnames(output) <- c("spread", "buy_signals", "sell_signals")
    return(output)
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
results <- calculate_beta_and_level(x, y, start_date, end_date)
results$beta
results$level
plot(results$spread, ylab = "Spread Value", main = "y - beta * x", cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

# Out of sample start and end dates
start_date_out_sample <- "2018-09-01"
end_date_out_sample <- "2019-09-01"
range <- paste(start_date_out_sample, "::", end_date_out_sample, sep = "")
# Out of sample analysis
spread_out_of_sample <- calculate_spread(x[range], y[range], results$beta)
plot(spread_out_of_sample, main = "y - beta * x", cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)
abline(h = results$level, lwd = 2)
