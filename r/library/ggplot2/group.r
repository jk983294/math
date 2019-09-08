library(ggplot2)
library(car)
ggplot(data = Salaries, aes(x = salary, fill = rank)) + geom_density(alpha = 0.3)

ggplot(Salaries, aes(x = yrs.since.phd, y = salary, color = rank, shape = sex)) + 
    geom_point()
