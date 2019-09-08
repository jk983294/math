m2 <- matrix(c(3:14), nrow = 4, byrow = TRUE, dimnames = list(c("row1", "row2", "row3", 
    "row4"), c("col1", "col2", "col3")))

# Accessing Elements of a Matrix
m2[1, 3]  # get 5
m2[2, ]  # access only the 2nd row, [6    7    8]
m2[, 3]  # access only the 3rd column, [5    8   11   14]^T
