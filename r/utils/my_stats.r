mystats <- function(x, na.omit = FALSE) {
    if (na.omit) 
        x <- x[!is.na(x)]
    m <- mean(x)
    n <- length(x)
    s <- sd(x)
    skew <- sum((x - m)^3/s^3)/n
    kurt <- sum((x - m)^4/s^4)/n - 3
    return(c(n = n, mean = m, stdev = s, skew = skew, kurtosis = kurt))
}

get_log_return <- function(prices) {
    diff(log(prices))
}

remove_missing_rows <- function(data) {
    # data is data frame or xts object
    return(data[complete.cases(data), ])
}
