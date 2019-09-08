m1 <- matrix(c(3:14), nrow = 4, byrow = TRUE)
m2 <- matrix(c(3:14), nrow = 4, byrow = TRUE, dimnames = list(c("row1", "row2", "row3", "row4"), c("col1", "col2", "col3")))

length(m1)  # 12 matrix is actually vector with dimension info
dim(m1)  # [1] 4 3
nrow(m1)  # 4
ncol(m1)  # 3

# equals
m1 == m2  # get element by element compare, return a matrix with bool elements
all(m1 == m2)  # TRUE
identical(m1, m2)  # FALSE

class(m1)  # 'matrix'
attributes(m1)  # $dim [1] 4 3
