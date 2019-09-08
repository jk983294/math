library(ggplot2)
data(singer, package = "lattice")
ggplot(singer, aes(x = height)) + geom_histogram()
