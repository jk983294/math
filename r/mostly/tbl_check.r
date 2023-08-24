library(datasets)
library(psych)

class(iris)
typeof(iris)

head(iris, 5)
tail(iris, 5)
View(iris) # invoke a spreadsheet-style data viewer

summary(iris)
psych::describe(iris)

(dt <- data.table(group = rep(letters[1:3], each = 5), value = 1:15))
dt[order(-value)][ , head(.SD, 3), by = group]  # max n data in each group
