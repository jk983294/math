library(ggplot2)
library(car)
ggplot(mtcars, aes(x = wt, y = mpg, size = disp)) + geom_point(shape = 21, color = "black", 
    fill = "cornsilk") + labs(x = "Weight", y = "Miles Per Gallon", title = "Bubble Chart", 
    size = "Engine\nDisplacement")
