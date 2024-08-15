x <- c(1., 2, 3, 4, 5, NA)
y <- c(6., 7, 8, 9, 10, NA)

# null hypothesis: the variances are equal
sd(x, na.rm = TRUE)
sd(y, na.rm = TRUE)
result <- var.test(x, y)

print(result)
