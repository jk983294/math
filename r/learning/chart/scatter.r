df <- data.frame(matrix(rnorm(20, 1), ncol = 2))
plot(x = df$X1, y = df$X2)
plot(x = df$X1, y = df$X2,
  pch = 19,         # Solid circle
  cex = 1.5,        # Make 150% size
  col = "#cc0000",  # Red
  main = "x vs y",
  xlab = "x",
  ylab = "y")

# find the correlation between one variable versus the remaining ones
df1 <- data.frame(matrix(rnorm(40, 1), ncol = 4))
pairs(~X1 + X2 + X3 + X4, data = df1, main = "Scatter Matrix")

library(car)
scatterplotMatrix(~mpg + disp + drat + wt, data = mtcars, spread = FALSE, smoother.args = list(lty = 2), 
    main = "Scatter Plot Matrix via car Package")

# 3D Scatterplots
library(scatterplot3d)
attach(mtcars)
scatterplot3d(wt, disp, mpg, main = "Basic 3D Scatter Plot")

scatterplot3d(wt, disp, mpg, pch = 16, highlight.3d = TRUE, type = "h", main = "3D Scatter Plot with Vertical Lines")

s3d <- scatterplot3d(wt, disp, mpg, pch = 16, highlight.3d = TRUE, type = "h", main = "3D Scatter Plot with Vertical Lines and Regression Plane")
fit <- lm(mpg ~ wt + disp)
s3d$plane3d(fit)
detach(mtcars)

# spinning 3D plot
library(rgl)
attach(mtcars)
plot3d(wt, disp, mpg, col = "red", size = 5)

library(car)
with(mtcars, scatter3d(wt, disp, mpg))
