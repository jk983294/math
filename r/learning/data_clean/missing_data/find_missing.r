# load the dataset
data(sleep, package = "VIM")


# list the rows that do not have missing values
sleep[complete.cases(sleep), ]


# list the rows that have one or more missing values
sleep[!complete.cases(sleep), ]


# tabulate missing values patterns
library(mice)
md.pattern(sleep)  # 第一列表示缺失个数,最后一列缺失变量个数


# plot missing values patterns
library("VIM")
aggr(sleep, prop = FALSE, numbers = TRUE)  # 每个变量的缺失值数,每个变量组合的缺失值数
matrixplot(sleep) # 每个实例数据的图形, 缺失值为红色
marginplot(sleep[c("Gest", "Dream")], pch = c(20), col = c("darkgray", "red", "blue"))


# use correlations to explore missing values
x <- as.data.frame(abs(is.na(sleep)))
head(sleep, n = 5)
head(x, n = 5)
y <- x[which(apply(x, 2, sum) > 0)]
cor(y) # Dream 和 NonD 常常一起缺失
cor(sleep, y, use = "pairwise.complete.obs")
