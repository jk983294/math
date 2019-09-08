mydata <- data.frame(matrix(rnorm(20, 1), ncol = 2))

# update whole column
mydata$sum1 <- mydata$X1 + mydata$X2  # method 1

attach(mydata)
mydata$sum2 <- X1 + X2  # method 2
detach(mydata)

mydata <- transform(mydata, sum3 = X1 + X2)  # method 3
transform(mydata, myvar = scale(X1) * 10 + 50)  # scale X1 by mean 50, sd = 10

# conditional update, variable[condition] <- expression
mydata <- within(mydata, {
    var <- NA
    var[X1 > 1] <- "big"
    var[X1 >= -1 & X1 <= 1] <- "middle"
    var[X1 < -1] <- "small"
})

mydata$x4 <- ifelse(mydata$X1 > 3, 1, 0)  # method 4
