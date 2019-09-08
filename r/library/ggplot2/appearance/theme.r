library(ggplot2)
library(car)
mytheme <- theme(plot.title = element_text(face = "bold.italic", size = "14", color = "brown"), 
    axis.title = element_text(face = "bold.italic", size = 10, color = "brown"), 
    axis.text = element_text(face = "bold", size = 9, color = "darkblue"), panel.background = element_rect(fill = "white", 
        color = "darkblue"), panel.grid.major.y = element_line(color = "grey", linetype = 1), 
    panel.grid.minor.y = element_line(color = "grey", linetype = 2), panel.grid.minor.x = element_blank(), 
    legend.position = "top")

ggplot(Salaries, aes(x = rank, y = salary, fill = sex)) + geom_boxplot() + labs(title = "Salary by Rank and Sex", 
    x = "Rank", y = "Salary") + mytheme
