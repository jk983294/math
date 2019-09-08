library(ggplot2)
data(singer, package = "lattice")
ggplot(singer, aes(x = voice.part, y = height)) + geom_violin(fill = "lightblue") + 
    geom_boxplot(fill = "lightgreen", width = 0.2)
