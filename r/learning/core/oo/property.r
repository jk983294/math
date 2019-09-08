# name property
die <- 1:6
names(die) <- c("one", "two", "three", "four", "five", "six")
names(die)
attributes(die) # $names

# dimension property
dim(die) <- c(2, 3)
dim(die) <- c(3, 2)

length(die)  # 6
dim(die)  # (3, 2)
str(die)  # 显示某个对象的结构 int [1:3, 1:2] 1 2 3 4 5 6
class(die)  # 显示某个对象的类 "matrix"
mode(die)  # 显示某个对象的模式 numeric"
names(die)  # 显示某对象中各成分的名称 NULL
edit(matrix) # investigate matrix code
