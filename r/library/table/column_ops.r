library("data.table")

(data <- data.table(x1 = 1:5, x2 = 6:10, x3 = 11:15))
(dt <- data.table(x1 = 1:5, x2 = letters[6:2], x3 = 9:5))

# display all column names
(cns <- colnames(dt))

# select and rename column
(ans <- data[, .(x4 = x1, x5 = x3)])

# in-place rename column
(setnames(data, c("x1", "x3"), c("x4", "x5"))) # x1 -> x4, x3 -> x5
(colnames(dt)[1L] <- "new_x1")

# in-place remove columns
data[, c("x2", "x5") := NULL]

# change column type
sapply(dt, class)
dt[, x1 := as.character(x1)] # one column
change_columns <- c("x1", "x3")
dt[, (change_columns) := lapply(.SD, as.character), .SDcols = change_columns]
dt[, (c("x1", "x3")) := lapply(.SD, as.character), .SDcols = c("x1", "x3")]
