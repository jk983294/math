library(ggplot2)

x <- 1:30
y <- c(5, 4, 5, 3, 4, 6, 7, 4, 2, 1, 5, 6, 2, 2, 2, 5, 6, 7, 8, 6, 2, 8, 7, 8, 8, 3, 3, 7, 2, 1)
z <- c("A", "B", "A", "A", "B", "A", "B", "B", "A", "A", "B", "A", "A", "A", "B", "A", "A", "B", "A", "B", "B", "A", "A", "A", "B", 
    "A", "A", "B", "A", "B")
df <- data.frame(x, y, z)

ggplot(df) + aes(x, y, color = z) + geom_point() + theme_grey()
ggplot(df) + aes(x, y, color = z) + geom_point() + theme_bw()
ggplot(df) + aes(x, y, color = z) + geom_point() + theme_minimal()
ggplot(df) + aes(x, y, color = z) + geom_point() + theme_classic()

# advanced theme
library(ggthemes)

ggplot(df) + aes(x, y, color = z) + geom_point() + theme_stata() + scale_colour_stata()
ggplot(df) + aes(x, y, color = z) + geom_point() + theme_excel() + scale_colour_excel()
ggplot(df) + aes(x, y, color = z) + geom_point() + theme_economist() + scale_colour_economist()
# fancy
ggplot(df) + geom_line(aes(x, y, colour = z), size = 1.5) + theme_economist() + scale_colour_economist() + theme(legend.position = "bottom", 
    axis.title = element_text(size = 12), legend.text = element_text(size = 9), legend.title = element_text(face = "bold", size = 9)) + 
    ggtitle("Title")
