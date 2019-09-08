library(ggplot2)

x <- 1:30
y <- c(5, 4, 5, 3, 4, 6, 7, 4, 2, 1, 5, 6, 2, 2, 2, 5, 6, 7, 8, 6, 2, 8, 7, 8, 8, 3, 3, 7, 2, 1)
z <- c("A", "B", "A", "A", "B", "A", "B", "B", "A", "A", "B", "A", "A", "A", "B", "A", "A", "B", "A", "B", "B", "A", "A", "A", "B", 
    "A", "A", "B", "A", "B")
df <- data.frame(x, y, z)

qplot(x, y, data = df, geom = "point")  # scatter
qplot(z, data = df, geom = "bar")
qplot(z, x, data = df, geom = "boxplot")
qplot(x, y, data = df, geom = "line")
qplot(y, data = df, geom = "histogram", binwidth = 3)
qplot(y, data = df, geom = "density")  # similar to histogram but function is smoothed

# corresponding ggplot
ggplot(df) + aes(x, y) + geom_point()
ggplot(df) + aes(z) + geom_bar()
ggplot(df) + aes(z, x) + geom_boxplot()
ggplot(df) + aes(x, y) + geom_line()
ggplot(df) + aes(y) + geom_histogram(binwidth = 3)
ggplot(df) + aes(y) + geom_density()
ggplot(df) + aes(y) + geom_freqpoly(binwidth = 1)
