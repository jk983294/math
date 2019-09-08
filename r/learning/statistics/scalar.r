x <- c(12, 7, 3, 4.2, 18, 2, 54, -21, 8, -5)
mean(x)  # 8.22
mean(x, trim = 0.3)  # 5.55 截尾平均数, drop 30% max and min values, then calculate the mean
mean(x, na.rm = TRUE)  # 8.22 drop NA values
sd(x)  # 19.20057 standard deviation
var(x)  # 368.6618 方差
median(x)  # 5.6
mad(x)  # 7.413 绝对中位差(median absolute deviation)
sum(x)  # 82.2
max(x)  # 54
min(x)  # -21
range(x)  # 求值域 [-21, 54]
IQR(x) # describes the middle 50% of values, difference between Q3 and Q1
quantile(x, 0.3)  # 2.7
quantile(x, c(0.3, 0.6))  # [2.7, 7.4]
diff(x)  # 差分 [-5.0  -4.0   1.2  13.8 -16.0  52.0 -75.0  29.0 -13.0]
diff(x, lag = 2)  # [-9.0  -2.8  15.0  -2.2  36.0 -23.0 -46.0  16.0]
scale(x)  # scale to Gaussian dist, mean = 0, sd = 1
scale(x) * 2 + 3  # scale to Gaussian dist, mean = 3, sd = 2
summary(x) # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.

# mode is value that has highest number of occurrences in a set of data
getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}
charv <- c("o", "it", "the", "it", "it")
getmode(charv)  # 'it'

# descriptive stats via sapply
mystats <- function(x, na.omit = FALSE) {
    if (na.omit) 
        x <- x[!is.na(x)]
    m <- mean(x)
    n <- length(x)
    s <- sd(x)
    skew <- sum((x - m)^3/s^3)/n
    kurt <- sum((x - m)^4/s^4)/n - 3
    return(c(n = n, mean = m, stdev = s, skew = skew, kurtosis = kurt))
}

df <- mtcars[c("mpg", "hp", "wt")]
sapply(df, mystats)

library(Hmisc)
Hmisc::describe(df)  # 观测的数量、缺失值和唯一值的数目、平均值、分位数,以及五个最大的值和五个最小的值


library(pastecs)
# 计算其中所有值、空值、缺失值的数量,以及最小值、最大值、值域,还有总和
# 中位数、平均数、平均数的标准误、平均数置信度为95%的置信区间、方差、标准差以及变异系数
# norm=TRUE则返回正态分布统计量,包括偏度和峰度(以及它们的统计显著程度)和Shapiro-Wilk正态检验结果
stat.desc(df)
