# cbind in R appends or combines vector, matrix or data frame by columns
df1 = data.frame(CustomerId = c(1:4), Product = c(rep("a", 2), rep("b", 2)))
df2 = data.frame(CustomerId = c(4:7), Product = c(rep("c", 2), rep("d", 2)))
rbind(df1, df2)

rbind(1, 1:3)  # short one will recycled
rbind(a = 1, b = 1:3)  # give name
