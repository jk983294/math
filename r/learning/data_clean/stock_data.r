source("utils/utils.r")

filter_valid_symbols <- function(symbols) {
    symbols <- toupper(symbols)
    valid <- regexpr("^[0-9]{6}.(SSE|SZE)$", symbols)
    # return only the valid ones
    return(sort(symbols[valid == 1]))
}

extract_prices <- function(filtered_symbols, file_path) {
    prices <- NULL
    for (i in 1:length(filtered_symbols)) {
        i_path <- paste(file_path, substr(filtered_symbols[i], 1, 6), ".csv", sep = "")
        dat <- read_stock_data(i_path)
        dat <- dat$close
        names(dat)[1] <- paste(substr(filtered_symbols[i], 1, 6), ".close", sep = "")
        if (i == 1) 
            prices <- dat else {
            prices <- merge(prices, dat)
        }
    }
    return(prices)
}

filter_prices <- function(prices) {
    # Returns a boolean vector of good or bad rows
    valid_rows <- complete.cases(prices)
    missing_rows <- which(valid_rows == FALSE)
    return(missing_rows)
}

compute_pairwise_correlations <- function(prices) {
    # calculates pairwise correlations of returns and plots the pairwise relationships
    returns <- apply(prices, 2, function(x) diff(log(x)))  # convert prices to returns
    pairs(returns, main = "Pairwise return scatter plot")
    correlation_matrix <- cor(returns, use = "complete.obs")
    return(correlation_matrix)
}

file_path <- "../data/"
symbols <- filter_valid_symbols(c("600848.SSE", "000001.SZE"))
prices <- extract_prices(symbols, file_path)
missing_rows <- filter_prices(prices)
prices <- prices[-missing_rows]
correlation_matrix <- compute_pairwise_correlations(prices)

## handle data
class(stk)  # read from csv, it is data.frame
stk <- aapl[rev(rownames(stk)), , drop = FALSE]  # reverse whole df
