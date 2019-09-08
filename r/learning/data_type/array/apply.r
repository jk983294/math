array1 <- array(1:12, dim = c(2, 3, 2))

# Calculations Across Array Elements, apply(array, margin, function ), margin
# rule: 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns
apply(array1, c(1), sum)
