library("VIM")

# method 1
newdata1 <- na.omit(sleep)

# method 2
newdata2 <- sleep[complete.cases(sleep), ]

nrow(sleep) # 62
nrow(newdata1) # 42
nrow(newdata2) # 42
