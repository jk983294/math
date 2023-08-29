library(data.table)
library(moments)
library(dplyr)

dt <- data.table(x = c("a", "b", "a", "b"), y = c(2., 1., 3., NA), z = c(1., 3., NA, 2.), m = c("m", "m", "n", "m"))

# count
(ans <- dt[, .N]) # how many rows in dt
(ans <- dt[, sum((y + z) < 5)]) # how many rows which (y + z) < 5
(ans <- dt[x == "b", length(x)]) # how many rows which x == 'b'
(ans <- dt[x == "b", .N]) # how many rows which x == 'b' using special symbol .N
(ans <- dt[, .(.N), by = .(x)]) # each group count
(ans <- dt[, .N, x]) # count by group of x
dt[, .N, .(x, m)][order(x, -N)]

# mean of selected rows
(ans <- dt[x == "a", .(m_y = mean(y, na.rm = TRUE), m_z = mean(z))])

# group
(ans <- dt[, .(.N), by = .(x)]) # each group count
(ans <- dt[y > 1, .(.N), by = x]) # each group count of x which row y > 1
(ans <- dt[, .(.N), by = .(x, y)]) # each group count
(ans <- dt[, .(mean(y, na.rm = TRUE), mean(z)), by = .(x)]) # each group mean

# group min/max
dt[, .SD[which.max(y)], by = x]
dt[, .SD[which.min(y)], by = x]
dt[, .(max_y = max(y, na.rm = TRUE)), by = x]
dt[, .(min_y = min(y, na.rm = TRUE)), by = x]

# top n of group
dt[order(-y)][, head(.SD, 1), by = x] # top n value by x

# apply to all subset using .SD
dt[, print(.SD), by = x] # .SD contains all the columns except the grouping columns by default
dt[, head(.SD, 1), by = x] # first row for each x
# group function apply to all columns except by columns
dt[, lapply(.SD, mean, na.rm = TRUE), by = .(x, m)]

dt_stats <- function(x) {
  mm <- colMeans(x, na.rm = TRUE)
  ss <- sapply(x, sd, na.rm = TRUE)
  median_ <- sapply(x, median, na.rm = TRUE)
  min_ <- sapply(x, min, na.rm = TRUE)
  max_ <- sapply(x, max, na.rm = TRUE)
  skew_ <- sapply(x, skewness, na.rm = TRUE)
  kurtosis_ <- sapply(x, kurtosis, na.rm = TRUE)
  list(names = names(x), median = median_, mean = mm, sd = ss, min = min_, max = max_, skew = skew_, kurtosis = kurtosis_)
}

dt[, dt_stats(.SD), by = .(x, m)] # group stats
dt[, dt_stats(.SD), .SDcols = c("y", "z")] # whole column stats

# expression in by
(ans <- dt[, .(.N), by = .(y > 1.5)]) # each group count
(ans <- dt[, .(.N), by = .(x, condtion_y = y > 1.5)]) # each group count

# group and sort
# it retaining the original order of groups is intentional and by design.
# keyby means to automatically sort by the variables in our grouping.
# keyby is typically faster than by because it doesn’t require recovering the original table’s order step.
(ans <- dt[, .(mean(y), mean(z)), keyby = .(x)]) # each group mean and sort by key
(ans <- dt[, .(mean(y), mean(z)), by = .(x)][order(x)]) # each group mean and sort by key
(ans <- dt[, .(mean(y), mean(z)), by = .(x)][order(-x)]) # reverse sort