fit2 <- lm(weight ~ height + I(height^2), data = women)
summary(fit2)
plot(women$height, women$weight, main = "Women Age 30-39", xlab = "Height (in inches)", 
    ylab = "Weight (in lbs)")
lines(women$height, fitted(fit2))
