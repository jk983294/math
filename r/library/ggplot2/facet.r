data(mtcars)
mtcars$am <- factor(mtcars$am, levels = c(0, 1), labels = c("Automatic", "Manual"))
mtcars$vs <- factor(mtcars$vs, levels = c(0, 1), labels = c("V-Engine", "Straight Engine"))
mtcars$cyl <- factor(mtcars$cyl)


library(ggplot2)
# am和vs是刻面变量, cyl是分组变量
ggplot(data = mtcars, aes(x = hp, y = mpg, shape = cyl, color = cyl)) + geom_point(size = 3) + 
    facet_grid(am ~ vs) + labs(title = "Automobile Data by Engine Type", x = "Horsepower", 
    y = "Miles Per Gallon")

# examples
data(singer, package = "lattice")
library(ggplot2)
ggplot(data = singer, aes(x = height)) + geom_histogram() + facet_wrap(~voice.part, 
    nrow = 4)

library(ggplot2)
library(car)
ggplot(Salaries, aes(x = yrs.since.phd, y = salary, color = rank, shape = rank)) + 
    geom_point() + facet_grid(. ~ sex)

data(singer, package = "lattice")
library(ggplot2)
ggplot(data = singer, aes(x = height, fill = voice.part)) + geom_density() + facet_grid(voice.part ~ 
    .)
