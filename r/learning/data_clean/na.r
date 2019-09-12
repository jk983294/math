y <- c(1, 2, 3, NA)
is.na(y)  # [FALSE FALSE FALSE TRUE]

# NA calculation
x <- c(1, 2, NA, 3)
z <- sum(x)  # NA
y <- sum(x, na.rm = TRUE)  # 6

# NULL calculation
sum(c(1, 2, NULL, 3))  # 6
sum(c(1, 2, NA, 3))  # NA
x <- NA
length(x)  # 1
isTRUE(x > 0) # FALSE
x <- NULL
length(x)  # 0
isTRUE(x > 0) # FALSE

# excluding missing values from analyses
mydata <- read.table(file = "data/to_process.txt", sep = ",", header = TRUE, quote = "\"")
sum(is.na(mydata))  # how many missing values
mydata <- na.omit(mydata)  # exclude missing value rows
