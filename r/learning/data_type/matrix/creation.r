# Matrices elements are arranged in a two-dimensional rectangular layout.
# Elements are the same atomic types.
m1 <- matrix(c(3:14), nrow = 4, byrow = TRUE)
colnames(m1) <- c("V1", "V2", "V3")
m2 <- matrix(c(3:14), nrow = 4, byrow = TRUE, dimnames = list(c("row1", "row2", "row3", 
    "row4"), c("col1", "col2", "col3")))
