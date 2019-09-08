# 获取来自给定均值向量和协方差阵的多元正态分布的数据
library(MASS)
mean <- c(230.7, 146.7, 3.6) # 均值向量
sigma <- matrix(c(15360.8, 6721.2, -47.1, 6721.2, 4700.9, -16.5, -47.1, -16.5, 0.3), 
    nrow = 3, ncol = 3) # 协方差阵
mydata <- mvrnorm(500, mean, sigma)
mydata <- as.data.frame(mydata)
names(mydata) <- c("y", "x1", "x2")
dim(mydata)
head(mydata, n = 10)
