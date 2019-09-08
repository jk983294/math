library(vioplot)

# 小提琴图基本上是核密度图以镜像方式在箱线图上的叠加
# 在图中,白点是中位数,黑色盒型的范围是下四分位点到上四分位点,细黑线表示须。外部形状即为核密度估计
x1 <- mtcars$mpg[mtcars$cyl == 4]
x2 <- mtcars$mpg[mtcars$cyl == 6]
x3 <- mtcars$mpg[mtcars$cyl == 8]
vioplot(x1, x2, x3, names = c("4 cyl", "6 cyl", "8 cyl"), col = "gold")
title("Violin Plots of Miles Per Gallon")
