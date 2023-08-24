(dt <- data.table(x1 = 1:6, x2 = 7:2, x3 = rep(1:2, each = 3)))

# cumsum
dt[, cs_x1 := cumsum(x1), by = x3]

# shift, lag one column
dt[, lag1 := shift(x2, n = 1, fill = x2[1L]), by = x3] # fill using first
dt[, lag3 := shift(x2, n = 1, fill = x2[.N]), by = x3] # fill using last
dt[, lag1_m_x1 := x1 * shift(x2), by = x3]
dt[, lag3_m_x1 := x1 * shift(x2, n = 3), by = x3]

# lag multi columns
cols <- c("x1", "x2")
(anscols <- paste("lag", cols, sep = "_"))
dt[, (anscols) := shift(.SD, n = 1, fill = 0, type = "lag"), .SDcols = cols]

# create id number by group, useful for some str/factor column
dt[, id_num := .GRP, by = x3] # id_num := as.numeric(factor(dt[, x3]))