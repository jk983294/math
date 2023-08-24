dt <- data.table(x = c("a", "b", "c", "b"), y = c(2, 1, 3, 4), z = c(1, 3, 2, 2))

# in-place remove columns
(dt[, c("x2", "x3") := NULL])