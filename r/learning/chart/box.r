# it are a measure of how well distributed is the data in a data set.

# 最小值、下四分位数、中位数、上四分位数以及最大值,描述了连续型变量的分布

# basic, it represents group X1 by X2
m1 <- data.frame(cbind(1:12, rep(1:3, each = 4)))
boxplot(X1 ~ X2, data = m1)
boxplot(X1 ~ X2, data = m1, notch = TRUE)

# below means 1:12 belong to same group
boxplot(1:12)

# Box plots for two crossed factors
mtcars$cyl.f <- factor(mtcars$cyl, levels = c(4, 6, 8), labels = c("4", "6", "8"))  # create a factor for number of cylinders
mtcars$am.f <- factor(mtcars$am, levels = c(0, 1), labels = c("auto", "standard"))  # create a factor for transmission type

# generate boxplot, group by <am.f, cyl.f>
boxplot(mpg ~ am.f * cyl.f, data = mtcars, varwidth = TRUE, col = c("gold", "darkgreen"), 
    main = "MPG Distribution by Auto Type", xlab = "Auto Type")
