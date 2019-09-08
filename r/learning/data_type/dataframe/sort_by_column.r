df <- data.frame(matrix(rnorm(20, 1), ncol = 4))

df[order(df$X1), ] # ascending order
df[order(df$X1, -df$X2), ] # X1 ascending, X2 descending
