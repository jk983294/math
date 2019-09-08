df <- data.frame(matrix(rnorm(20, 1), ncol = 2))
df[sample(1:nrow(df), 5, replace = FALSE), ]  # sample 5 observations
