library(data.table)
library(moments)

dt <- data.table(x = c("a", "b", "a", "b"), y = c(2., 1., 3., NA), z = c(1., 3., NA, 2.))

ans <- dt[, sum((y + z) < 5)]  # how many rows which (y + z) < 5
ans

ans <- dt[x == 'b', length(x)]  # how many rows which x == 'b'
ans

ans <- dt[x == 'b', .N]  # how many rows which x == 'b' using special symbol .N
ans

# mean of selected rows
ans <- dt[x == 'a', .(m_y = mean(y), m_z = mean(z))]
ans

# group
ans <- dt[, .(.N), by = .(x)]   # each group count
ans
ans <- dt[y > 1, .(.N), by = x]   # each group count
ans
ans <- dt[, .(.N), by = .(x, y)]   # each group count
ans
ans <- dt[, .(mean(y), mean(z)), by = .(x)]   # each group mean
ans

# apply to all subset using .SD
dt[, print(.SD), by = x]  # .SD contains all the columns except the grouping columns by default
dt[, head(.SD, 1), by = x]  # first row for each x
# group function apply to all columns
dt[, lapply(.SD, mean), by = x]

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

dt[, dt_stats(.SD), by = x]                 # group stats
dt[, dt_stats(.SD), .SDcols = c("y", "z")]  # whole column stats

# expression in by
ans <- dt[, .(.N), by = .(y > 1.5)]   # each group count
ans

# group and sort
ans <- dt[, .(mean(y), mean(z)), keyby = .(x)]   # each group mean and sort by key
ans
ans <- dt[, .(mean(y), mean(z)), by = .(x)][order(x)]   # each group mean and sort by key
ans
ans <- dt[, .(mean(y), mean(z)), by = .(x)][order(-x)]   # reverse sort
ans
