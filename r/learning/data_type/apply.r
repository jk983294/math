df <- data.frame(matrix(rnorm(30, 1), ncol = 5))

# rapply - apply a function to all elements of a list recursively
rapply(df, function(x) {
    x + 1
}, class = c("numeric"))

# apply by column or row
apply(df, 1, sum)  # row wise sum up
apply(df, 2, sum)  # column wise sum up

# lapply function returns list as output
lapply(df, sum)  # column wise

# sapply similar to lapply function but returns vector as output
sapply(df, sum)  # column wise

# mapply is a multivariate version of sapply
mapply(sum, df)  # the same as sapply(df, sum)
mapply(sum, df, df)  # twice of sapply(df, sum)
