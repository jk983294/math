# think of .SD as a symbol representing “each group"
# .BY is a special symbol that holds the value of by
# .N is an integer containing the number of rows in the group.
# .I is an integer vector equal to seq_len(nrow(x))
# .GRP is an integer containing a simple group counter
# .NGRP is an integer containing the number of groups

help("special-symbols")

dt <- data.table(x = c("a", "b", "a", "b"), y = c(2., 1., 3., NA), z = c(1., 3., NA, 2.), m = c("m", "m", "n", "m"))

# print all column by group of x
dt[, print(.SD), by = x]

# find max y regard group of x and print all columns
dt[, .SD[which.max(y)], by = x]
dt[, .SD[which.min(y)], by = x]


(DT <- data.table(x = rep(c("b", "a", "c"), each = 3), v = c(1, 1, 1, 2, 2, 1, 1, 2, 2), y = c(1, 3, 6), a = 1:9, b = 9:1))
(X <- data.table(x = c("c", "b"), v = 8:7, foo = c(4, 2)))
DT[, .SD, .SDcols = x:y] # select columns between 'x' and 'y'
DT[, .SD[1]] # first row of all columns
DT[, .SD[1], by = x] # first row of 'y' and 'v' for each group in 'x'
DT[, c(.N, lapply(.SD, sum)), by = x] # get row number and sum of columns by group x
DT[, .I[1], by = x] # first row's number in DT by each group x
rleid(DT$v) # generate run-length type group id of 'v'
DT[, .N, by = rleid(v)] # get count of consecutive runs of 'v'
DT[, c(.(y = max(y)), lapply(.SD, min)),
    by = rleid(v), .SDcols = v:b
] # compute 'j' for each consecutive runs of 'v'
DT[, grp := .GRP, by = x] # add a group id by 'x'
DT[, grp_pct := .GRP / .NGRP, by = x] # add a group "progress" counter
X[, DT[.BY, y, on = "x"], by = x] # join within each group