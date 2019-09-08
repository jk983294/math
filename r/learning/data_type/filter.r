# return subsets of vectors, matrices or data frames which meet conditions
df <- data.frame(cbind(matrix(rnorm(30, 1), ncol = 5), c(1, seq(5))))
names(df)  # 'X1' 'X2' 'X3' 'X4' 'X5' 'X6'
subset(df, X6 > 2)  # selectrows which X6 column > 2
subset(df, X6 > 2 & X1 > 0.6)  # multi-condition
subset(df, X6 > 2, select = c(X1, X4))  # selected columns
subset(df, X6 > 2, select = X2:X5)  # selected columns

# which return position satisfy condition
x <- 1:5
which(x > 2)

first_position <- function(v, n) {
    return(which(v > n)[1])
}

first_position(x, 2)
