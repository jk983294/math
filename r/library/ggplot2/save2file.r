ggplot(data = mtcars, aes(x = mpg)) + geom_histogram()
ggsave(file = "mygraph.pdf")
