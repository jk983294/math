seq(1, 10, 2)  # [1 3 5 7 9]
rep(1:3, 2)  # [1 2 3 1 2 3]

# 创建美观的分割点
pretty(1:20, n = 2)  # 0 10 20
pretty(1:20, n = 10)  # [0  2  4  6  8 10 12 14 16 18 20]

# cut 将连续型变量分割为n个水平的因子
x <- rnorm(100)
c <- cut(x, breaks = -5:5)
summary(c)
table(c)
