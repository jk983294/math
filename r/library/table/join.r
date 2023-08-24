# J : it is a direct alias of list
# SJ : Sorted Join. The same as J() but additionally setkey() is called on all columns
# CJ : Cross Join.

(DT <- data.table(A = 5:1, B = letters[5:1]))
(setkey(DT, B)) # reorders table and marks it sorted
DT[J("b")] # returns the 2nd row
DT[list("b")] # same
DT[.("b")] # same using the dot alias for list

# CJ usage examples
(dt <- CJ(c(5, NA, 1), c(1, 3, 2))) # sorted and keyed data.table
(dt <- do.call(CJ, list(c(5, NA, 1), c(1, 3, 2)))) # same as above
(dt <- CJ(c(5, NA, 1), c(1, 3, 2), sorted = FALSE)) # same order as input, unkeyed

# CJ 'unique=' argument
x <- c(1, 1, 2)
y <- c(4, 6, 4)
CJ(x, y) # output columns are automatically named 'x' and 'y'
CJ(x, y, unique = TRUE) # unique(x) and unique(y) are computed automatically