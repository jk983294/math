# DT[, c("colA", "colB", ...) := list(valA, valB, ...)]
# DT[, colA := valA]
# DT[, `:=`(colA = valA, colB = valB, )]

dt <- data.table(x = c("a", "b", "a", "b"), y = c(2., 1., 3., NA), z = c(1., 3., NA, 2.), m = c("m", "m", "n", "m"))

# add column
dt[, y_z := y + 2. * z]
dt[, max_y_by_x_m := max(y), by = .(x, m)]
in_cols  <- c("y", "z")
out_cols <- c("out_y", "out_z")
dt[, c(out_cols) := lapply(.SD, max), by = .(x, m), .SDcols = in_cols]

# update existing column
# sub-assign by reference
dt[x == "a", y_z := y_z + 2.]  # select rows to update

# delete column
dt[, y_z := NULL]
dt[, c("max_y_by_x_m") := NULL]
dt[, c(out_cols) := NULL]

# convert
cols <- c("x", "m")
dt[, lapply(.SD, is.character), .SDcols = cols] # confirm that they're stored as character
lapply(dt[, ..cols], is.character) # confirm that they're stored as character
dt[, (cols) := lapply(.SD, factor), .SDcols = cols] # convert columns to factor
unique(dt[["x"]])

# copy
dt_deep <- copy(dt)
