library(car)
library(ggplot2)

# position effect
ggplot(Salaries, aes(x = rank, fill = sex)) + geom_bar(position = "stack") + labs(title = "position=\"stack\"")
ggplot(Salaries, aes(x = rank, fill = sex)) + geom_bar(position = "dodge") + labs(title = "position=\"dodge\"")
ggplot(Salaries, aes(x = rank, fill = sex)) + geom_bar(position = "fill") + labs(title = "position=\"fill\"")

# Placing options 通常来说,变量应该设在 aes() 函数内,分配常数应该在 aes() 函数外
ggplot(Salaries, aes(x = rank, fill = sex)) + geom_bar()
ggplot(Salaries, aes(x = rank)) + geom_bar(fill = "red")
ggplot(Salaries, aes(x = rank, fill = "red")) + geom_bar()
