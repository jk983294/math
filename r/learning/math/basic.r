seq(1, 5)  # [1 2 3 4 5]
mean(1:10)  # 5.5
sum(1:10)  # 55
cumsum(1:10)  # [1 3 6 10 15 21 28 36 45 55]
prod(1:5)  # 120
cumprod(1:5)  # [1 2 6 24 120]
max(1:10)  # 10
min(1:10)  # 1
which.max(1:10)  # position 10
which.min(1:10)  # position 1
diff(seq(1, 10, 2))  # [2 2 2 2]

## remove first/last element
(1:5)[-1]  # [2 3 4 5]
(1:5)[-5]  # [1 2 3 4]
diffs <- (1:5)[-1] - (1:5)[-5]  # [1 1 1 1]

# column wise min max
x <- c(3, 26, 122, 6)
y <- c(43, 2, 54, 8)
z <- c(9, 32, 1, 9)
# pair wise max/min
pmax(x, y, z)  # [43 32 122 9]
pmin(x, y, z)  # [3 2 1 6]

# work on scalar/vector
abs(-2)
sqrt(4)
ceiling(-3.475)  # -3
ceiling(3.475)  # 4
floor(-3.475)  # -4
floor(3.475)  # 3
trunc(-3.475)  # -3
trunc(3.475)  # -3
round(3.14)  # 3
round(-3.1415, digits = 2)  # -3.14 舍入为指定位的小数
signif(-3.1415, digits = 2)  # -3.1 舍入为指定的有效数字位数
exp(2) # e^2
log(10, base = 10)  # 1
log(2.71828)  # 1
log10(10)  # 1
factorial(3)  # 6

# also can apply on vector
abs(c(-2, 1, -3))  # [2 1 3]

# triangle function
cos(x)
sin(x)
tan(x)
acos(x)
asin(x)
atan(x)
cosh(x)
sinh(x)
tanh(x)
acosh(x)
asinh(x)
atanh(x)
