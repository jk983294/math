library(data.table)

dt <- data.table(x = c("a", "b", "c", "b"), y = c(2, 1, 3, 4), z = c(1, 3, 2, 2))

#       DT[      i,                  j,           by]
# SQL:  where | order by   select | update  group by
# Take DT, subset/reorder rows using i, then calculate j, grouped by by

# subset rows in i
(ans <- dt[x == "b" & y > 3])
(ans <- dt[y %between% c(1.5, 3.5)])
(ans <- dt[1L:2L]) # get the first two rows
dt[.N] # last row of all columns
dt[, .SD[1L]] # first row of all columns
dt[order(-y)][, head(.SD, 2L), by = x] # top n value by x

# subset rows in vec
vec_y_to_select <- c(1, 2)
dt[y %in% vec_y_to_select, ]
dt[y %chin% vec_y_to_select, ] # %chin% more optimized than original %in%
dt[is.element(y, vec_y_to_select), ]
dplyr::filter(dt, y %in% vec_y_to_select)
setDT(dt, key = "y")[J(vec_y_to_select)]


# subset columns in j
(ans <- dt$z) # return vector z in data.frame way
(ans <- dt[, z]) # return vector z
(ans <- dt[, list(z)]) # return table contains column z
(ans <- dt[, list(y, z)]) # return table contains column y&z
(ans <- dt[, .(y, z)]) # .() list sugar, return table
(ans <- dt[, c("y", "z")]) # refer to columns by string names
dt[, .SD, .SDcols = x:z] # select columns between 'x' and 'y'
dt[, .SD, .SDcols = c("x", "y")] # select columns

# Why two dots? Think of it like the two dots in a Unix command-line terminal that move you up one directory.
# Here, you’re moving up one namespace, from the environment inside data.table brackets up to the global environment.
select_cols <- c("y", "z")
(ans <- dt[, ..select_cols])

# with = FALSE disables the ability to refer to columns as if they are variables
# thereby switch to the 'data.frame mode' instead of 'data.table mode'
(ans <- dt[, select_cols, with = FALSE])

# select and rename column
(ans <- dt[, .(new_y = y, new_z = z)])

# column range select
ans <- dt[, y:z]

# opposite select column
(ans <- dt[, !c("y", "z")])
(ans <- dt[, -c("y", "z")])
(ans <- dt[, !(y:z)])
(ans <- dt[, -(y:z)])

# in-place remove columns
(dt[ , c("y", "z") := NULL])
