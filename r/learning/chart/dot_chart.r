dotchart(mtcars$mpg, labels = row.names(mtcars), cex = 0.7, main = "Gas Mileage for Car Models", 
    xlab = "Miles Per Gallon")


# Dot plot grouped, sorted, and colored
x <- mtcars[order(mtcars$mpg), ]
x$cyl <- factor(x$cyl)
x$color[x$cyl == 4] <- "red"
x$color[x$cyl == 6] <- "blue"
x$color[x$cyl == 8] <- "darkgreen"
dotchart(x$mpg, labels = row.names(x), cex = 0.7, pch = 19, groups = x$cyl, gcolor = "black", 
    color = x$color, main = "Gas Mileage for Car Models\ngrouped by cylinder", xlab = "Miles Per Gallon")
