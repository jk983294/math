# refer to plot.r

# another example
with(mtcars, {
    plot(wt, mpg)
    abline(lm(mpg ~ wt))
    title("Regression of MPG on Weight")
})

# type= options in the plot() and lines() functions
x <- c(1:5)
y <- c(1:5)
par(mfrow = c(2, 4))
types <- c("p", "l", "o", "b", "c", "s", "S", "h")
for (i in types) {
    plottitle <- paste("type=", i)
    plot(x, y, type = i, col = "red", lwd = 2, cex = 1, main = plottitle)
}
par(mfrow = c(1, 1)) # restore back parameter

# Line chart displaying the growth of 5 Orange trees over time
summary(Orange)
head(Orange)
Orange$Tree <- as.numeric(Orange$Tree)
ntrees <- max(Orange$Tree)
xrange <- range(Orange$age)
yrange <- range(Orange$circumference)
plot(xrange, yrange, type = "n", xlab = "Age (days)", ylab = "Circumference (mm)")
colors <- rainbow(ntrees)
linetype <- c(1:ntrees)
plotchar <- seq(18, 18 + ntrees, 1)
for (i in 1:ntrees) {
    tree <- subset(Orange, Tree == i)
    lines(tree$age, tree$circumference,
        type = "b", lwd = 2, lty = linetype[i], col = colors[i],
        pch = plotchar[i]
    )
}
title("Tree Growth", "example of line plot")
legend(xrange[1], yrange[2], 1:ntrees,
    cex = 0.8, col = colors, pch = plotchar, lty = linetype,
    title = "Tree"
)