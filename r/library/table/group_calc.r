DT <- data.table(
  ID = c("b", "b", "b", "a", "a", "c"),
  a = 1L:6L,
  b = 7L:12L,
  c = 13L:18L
)

c(DT$a, DT$b)
DT[, .(val = c(a, b)), by = ID]
DT[, .(val = list(c(a, b))), by = ID] # for each group, we return a list of all concatenated values

DT[, print(c(a, b)), by = ID] # a vector is returned, with length = 6,4,2
DT[, print(list(c(a, b))), by = ID] # returns a list of length 1 for each group
