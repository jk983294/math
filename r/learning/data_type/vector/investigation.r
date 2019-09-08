# vector, hold elements of different classes, c stands for concatenation
v <- c(1, 2, 3)
is.vector(v)  # TRUE
class(v)  # 'numeric'
typeof(v)  # 'double'
str(v)  # num [1:3] 1 2 3
summary(v)  # min/max/mean

# property
die <- 1:6
is.vector(die)  # TRUE
length(die)  # 6
typeof(die)  # 'integer'
names(die)  # NULL
attributes(die)  # NULL

# check if element in vector
print(8 %in% 1:10)
print(18 %in% 1:10)
