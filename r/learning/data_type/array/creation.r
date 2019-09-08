# Arrays can store data in more than two dimensions.  dimension (2, 3, 4) means
# it creates 4 rectangular matrices each with 2 rows and 3 columns

## if element < total count, then recycling reused
array(c(5, 9, 3, 10, 11, 12, 13, 14, 15), dim = c(3, 3, 2))
a <- array(data = 1:8, dim = c(2, 2, 2))
a <- array(c("green", "yellow"), dim = c(3, 3, 2))

vector1 <- 1:5
vector2 <- 2:14
column.names <- c("COL1", "COL2", "COL3")
row.names <- c("ROW1", "ROW2", "ROW3")
matrix.names <- c("Matrix1", "Matrix2")
result <- array(c(vector1, vector2), dim = c(3, 3, 2), dimnames = list(column.names, 
    row.names, matrix.names))
result
