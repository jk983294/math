(DT1 <- data.table(B = letters[5L:1L], A = 5L:1L, C = 15L:11L))
(DT2 <- data.table(B = letters[3L:7L], C = 11L:15L, D = 21L:25L))
(setkey(DT1, B))
(setkey(DT2, B))

#### outer join
# take the rows of DT2 and add the information of DT1 for these rows
# keep all DT2's key column or B column
DT1[DT2]
DT1[DT2, on = c("B")]
merge(DT1, DT2, all.y = TRUE) # LEFT OUTER
DT2[, DT1[.SD, A, on = c("B")]] # DT2's key join with DT1 to get A column

# keep all DT1's key column or B column
DT2[DT1]
DT2[DT1, on = c("B")]
merge(DT1, DT2, all.x = TRUE) # RIGHT OUTER
DT1[, DT2[.SD, C, on = c("B")]] # DT1's key join with DT2 to get C column

#### inner join
DT1[DT2, on = c("B"), nomatch = 0]
DT2[DT1, on = c("B"), nomatch = 0]
merge(DT1, DT2, all = FALSE)

#### full outer Join
merge(DT1, DT2, all = TRUE)


# SJ : Sorted Join. The same as J() but additionally setkey() is called on all columns
# CJ : Cross Join.
(dt <- CJ(c(5L, NA, 1L), c(1L, 3L, 2L))) # sorted and keyed data.table
(dt <- do.call(CJ, list(c(5L, NA, 1L), c(1L, 3L, 2L)))) # same as above
(dt <- CJ(c(5L, NA, 1L), c(1L, 3L, 2L), sorted = FALSE)) # same order as input, unkeyed
x <- c(1L, 1L, 2L)
y <- c(4L, 6L, 4L)
CJ(x, y) # output columns are automatically named 'x' and 'y'
CJ(x, y, unique = TRUE) # unique(x) and unique(y) are computed automatically