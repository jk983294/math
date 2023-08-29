library(data.table)

dt <- data.table(school = c("NTU", "SMU", "NUS"), rank = c(2, 1, 3), size = c(1, 3, 2))
dim(dt)
class(dt$school) # character column never converted to factors by default

# add row
dt_row <- rbind(dt, list("a", 5, 6))
new_row <- data.table("school" = "a", "rank" = 5, size = 6)
dt_add_row <- rbindlist(list(dt, new_row))

# add column
new_column_name <- "ncn"
dt[, (new_column_name) := rank + 5]
dt[, "x_new" := NA]
dt[, `:=`(cast_demo, as.character(size))]
dt[, ("cast_demo1") := as.character(size)]
dt <- tibble::rownames_to_column(dt, "row_names") # add rownames column
dt$r_name <- row.names(dt) # add rownames column
# cbind, merge two column
cbind(dt, merged_column = paste0(dt[, size], c("-"), dt[, school]))

# add table
dt1 <- data.table(field1 = c(1, 1, 1), field2 = c(2, 2, 2), field3 = c(3, 3, 3))
dt_tbl <- cbind(dt, dt1[, 1:2])

# empty dt
data <- data.table(va = numeric(), vb = numeric(), vc = numeric())

# create from list
l <- list(1L:3L, 4L:6L, 7L:9L)
dt <- setDT(lapply(l, unlist))

# SJ : Sorted Join. The same as J() but additionally setkey() is called on all columns
# CJ : Cross Join.
(dt <- CJ(c(5L, NA, 1L), c(1L, 3L, 2L))) # sorted and keyed data.table
(dt <- do.call(CJ, list(c(5L, NA, 1L), c(1L, 3L, 2L)))) # same as above
(dt <- CJ(c(5L, NA, 1L), c(1L, 3L, 2L), sorted = FALSE)) # same order as input, unkeyed
x <- c(1L, 1L, 2L)
y <- c(4L, 6L, 4L)
CJ(x, y) # output columns are automatically named 'x' and 'y'
CJ(x, y, unique = TRUE) # unique(x) and unique(y) are computed automatically
