library(dplyr)

df <- data.frame(ID = c("a", "b", "a", "b", "a", "b"), x = c(1, 2, 3, 4, 5, 6), y = c(1, 4, 2, 3, 6, 2))
newdf <- group_by(df, ID)
summarize(newdf, count = n(), xbar = mean(x), yvar = var(y))
