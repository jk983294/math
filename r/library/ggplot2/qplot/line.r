library(ggplot2)
library(data.table)
library(dplyr)

# dot graph
x <- c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)
y <- x^3
y1 <- x^2
dt <- data.table(x, y, y1)

qplot(x, y)
qplot(x, y, geom = "line")

# multiple line method1
ggplot(dt, aes(x)) +
    geom_line(aes(y = y, colour = "y")) +
    geom_line(aes(y = y1, colour = "y1"))

# multiple line method2
melted = melt(dt, id.vars="x")
ggplot(data=melted, aes(x=x, y=value, group=variable, colour=variable)) + geom_line()
