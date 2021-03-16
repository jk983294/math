library(ggplot2)

data(singer, package = "lattice")
ggplot(singer, aes(x = height)) + geom_histogram()

scores <- 1:10
qplot(scores, geom = "histogram")
