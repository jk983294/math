x <- c(21, 62, 10, 53)
labels <- c("London", "New York", "Singapore", "Mumbai")

# naive pie chart
pie(x, labels)

# pie with percentage
piepercent <- round(100 * x/sum(x), 1)
pie(x, labels = piepercent, main = "City pie chart", col = rainbow(length(x)))
legend("topright", labels, cex = 0.8, fill = rainbow(length(x)))

# 3d pie charts
library(plotrix)
pie3D(x, labels = labels, explode = 0.1, main = "Pie Chart of Countries")

# 扇形图
library(plotrix)
fan.plot(x, labels = labels, main = "Fan Plot")
