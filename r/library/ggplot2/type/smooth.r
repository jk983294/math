library(ggplot2)
library(car)
ggplot(data = Salaries, aes(x = yrs.since.phd, y = salary)) + geom_smooth() + geom_point()

ggplot(data = Salaries, aes(x = yrs.since.phd, y = salary, linetype = sex, shape = sex, 
    color = sex)) + geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE, 
    size = 1) + geom_point(size = 2)
