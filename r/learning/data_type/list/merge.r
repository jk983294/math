list1 <- list(1:2)
list2 <- list(10:11)
merged1 <- c(list1, list2)  # concate, [[1, 2], [10, 11]]
sliced <- merged1[2]  # sliced = list2
merged2 <- merge(list1, list2)  # cartesian product
