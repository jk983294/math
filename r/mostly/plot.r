library(data.table)
library(ggplot2)

v <- runif(10)
plot(v, type = 'l')
ts.plot(v)

dt <- as.data.table(list(x=1L:10L, b=v, c=runif(10)))

ts.plot(cbind(dt$b, dt$c), gpars= list(col=rainbow(2)))

ggplot(dt, aes(x = x)) +
    geom_line(aes(y = b), color = "red") +
    geom_line(aes(y = c), color = "blue")
