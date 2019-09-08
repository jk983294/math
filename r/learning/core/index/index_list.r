set.seed(1234)
fit <- kmeans(iris[1:4], 3)
names(fit)
unclass(fit)
sapply(fit, class)

# Indexing lists
fit[c(2, 3)] # 成分in list
fit[2] # 成分in list
fit[[2]] # 得到成分中的元素
fit$centers # 命名成分, == fit[[2]]
fit[[2]][1, ] # first row of fit[[2]]
fit$centers$Petal.Width  # should give an error, $ operator is invalid for atomic vectors
