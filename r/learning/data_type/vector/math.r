# dot product
die <- 1:2  # [1, 2]
die %*% die  # 1*1 5

# cross product
die %o% die  # 2*2 [[1,2],[2,4]]

# element recycling
v1 <- c(3, 8, 4, 5, 0, 11)
v2 <- c(4, 11)
v1 + v2  # the shorter vector are recycled to complete the operations

# check condition
all(v1 > 10)  # FALSE
any(v1 > 10)  # TRUE

lag(v1)  # NA  3  8  4  5  0
lead(v1)  # 8  4  5  0 11 NA
cumsum(v1) # 3 11 15 20 20 31
cummean(v1) # 3.000 5.500 5.000 5.000 4.000 5.167
min_rank(v1) # 2 5 3 4 1 6
min_rank(desc(v1)) # 5 2 4 3 6 1
row_number(v1) # 2 5 3 4 1 6
dense_rank(v1) # 2 5 3 4 1 6
percent_rank(v1) # 0.2 0.8 0.4 0.6 0.0 1.0
cume_dist(v1) # 0.3333 0.8333 0.5000 0.6667 0.1667 1.0000
