library(ggplot2)
library(car)
ggplot(data = Salaries, aes(x = yrs.since.phd, y = salary, color = rank)) + scale_color_manual(values = c("orange", 
    "olivedrab", "navy")) + geom_point(size = 2)

ggplot(data = Salaries, aes(x = yrs.since.phd, y = salary, color = rank)) + scale_color_brewer(palette = "Set1") + 
    geom_point(size = 2)

library(RColorBrewer)
display.brewer.all()
