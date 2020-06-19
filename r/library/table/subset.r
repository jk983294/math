library(data.table)

dt <- data.table(x = c("a", "b", "c", "b"), y = c(2, 1, 3, 4), z = c(1, 3, 2, 2))

#       DT[      i,                  j,           by]
# SQL:  where | order by   select | update  group by
# Take DT, subset/reorder rows using i, then calculate j, grouped by by

# subset rows in i
ans <- dt[x == "b" & y > 3]
ans
ans <- dt[1:2]   # get the first two rows
ans

# subset columns in j
ans <- dt[, z]          # return vector
ans
ans <- dt[, list(z)]    # return table
ans
ans <- dt[, list(y, z)]
ans
ans <- dt[, .(y, z)]    # list sugar
ans
ans <- dt[, c('y', 'z')]    # refer to columns by string names 
ans
select_cols = c("y", "z")
ans <- dt[, ..select_cols]    # string name sugar
ans
# with = FALSE disables the ability to refer to columns as if they are variables, thereby restoring the “data.frame mode”
ans <- dt[, select_cols, with = FALSE]
ans
ans <- dt[, .(new_y = y, new_z = z)]    # rename column
ans

# column range select
ans <- dt[, y:z]

# opposite select column
ans <- dt[, !c('y', 'z')]
ans <- dt[, -c('y', 'z')]
ans <- dt[, !(y:z)]
ans <- dt[, -(y:z)]
