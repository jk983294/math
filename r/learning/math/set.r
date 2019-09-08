x <- c(1, 2, 5)
y <- c(5, 1, 8, 9)
union(x, y) # [1 2 5 8 9]
intersect(x, y) # [1 5]
setdiff(x, y) # [2]
setdiff(y, x) # [8 9]
setequal(x, y) # FALSE
setequal(x, c(1, 2, 5)) # TRUE
2 %in% x # TRUE
2 %in% y # FALSE

choose(3, 2) # 3
combn(1:3, 2) # column wise output, 3 columns
combn(1:3, 2, sum) # [3 4 5] apply sum to each result column
