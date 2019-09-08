# sorting
v <- c(3, 8, 4, 5, 0, 11, -9, 304)
sort(v) # [-9 0 3 4 5 8 11 304]
sort(v, decreasing = TRUE) # [304 11 8 5 4 3 0 -9]

# get original order
order(v) # [7 5 1 3 4 2 6 8]
order(v, decreasing = TRUE) # [8 6 2 4 3 1 5 7]
