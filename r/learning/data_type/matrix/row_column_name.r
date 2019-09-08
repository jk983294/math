m1 <- matrix(rnorm(16), ncol = 4)

dimnames(m1) <- list(rownames(rmatrix, do.NULL = FALSE, prefix = "row"), colnames(rmatrix, 
    do.NULL = FALSE, prefix = "col"))

m1

colnames(m1) <- c("C1", "C2", "C3", "C4")  # change column name
rownames(m1) <- c("R1", "R2", "R3", "R4")  # change row name
m1
