# mutate() adds new variables that are functions of existing variables
library(dplyr)

df <- data.frame(product = c("a", "b", "c"), revenue = c(10, 5, 7), cost = c(5, 4, 5))
mutate(df, profit = revenue - cost, profit_margin = profit/revenue)
mutate(df, a = ifelse(revenue > 6, 1, 0))
mutate(df, b = case_when(revenue == 10 | product == "b" ~ 0, revenue == 7 ~ 1))

# inplace mutate
transmute(df, a = ifelse(revenue > 6, 1, 0))
