# cbind in R appends or combines vector, matrix or data frame by columns
df1 = data.frame(x = c("a", "b", "c", "d"), y = c(1, 2, 3, 4))
df2 = data.frame(z = c("e", "f", "g", "h"), m = c(5, 6, 7, 8))
cbind(df1, df2)

m <- cbind(1:2, 3:4)
cbind(m, 8:9)[, c(1, 3, 2)]  # insert a column in between

cbind(0, 1:3)  # short one will recycled
cbind(a = 1, b = 1:3)  # give name

xx <- data.frame(I = rep(0, 2))
cbind(xx, X = rbind(a = 1, b = 1:3))  # named differently
