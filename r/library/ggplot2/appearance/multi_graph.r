library(ggplot2)
library(car)
p1 <- ggplot(data = Salaries, aes(x = rank)) + geom_bar()
p2 <- ggplot(data = Salaries, aes(x = sex)) + geom_bar()
p3 <- ggplot(data = Salaries, aes(x = yrs.since.phd, y = salary)) + geom_point()

library(gridExtra)
grid.arrange(p1, p2, p3, ncol = 3)
