dt <- as.data.table(matrix(rnorm(30L, 1L), ncol = 5L))

# lapply(): Loop over a list and evaluate a function on each element, return list
# sapply(): Same as lapply but try to simplify the result
# apply(): Apply a function over the margins of an array
# mapply(): Multivariate version of lapply

lapply(dt, sum)  # column wise, return list
sapply(dt, sum)  # column wise, return vec/matrix
apply(dt, 1L, sum)  # row wise sum up
apply(dt, 2L, sum)  # column wise sum up
apply(dt, 2L, quantile, probs = c(0.25, 0.75)) # column quantile
rowSums(dt) # = apply(x, 1, sum)
rowMeans(dt) # = apply(x, 1, mean)
colSums(dt) # = apply(x, 2, sum)
colMeans(dt) # = apply(x, 2, mean)

lapply(1L:3L, function(x, to_add) x + 1L + to_add, to_add = 3L)

mapply(rep, 1L:4L, 4L:1L) # = list(rep(1, 4), rep(2, 3), rep(3, 2), rep(4, 1))
mapply(sum, dt)  # the same as sapply(dt, sum)
mapply(sum, dt, dt)  # twice of sapply(dt, sum)

x <- c(1L:10L)
results <- numeric(length(x)) # pre-allocation
for (i in 1L:length(x)) results[i] <- log(x[i])
