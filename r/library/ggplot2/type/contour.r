library(ggplot2)
ggplot(faithfuld, aes(waiting, eruptions, z = density)) + geom_contour()
