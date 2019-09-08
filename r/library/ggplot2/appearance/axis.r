library(ggplot2)
library(car)
ggplot(data = Salaries, aes(x = rank, y = salary, fill = sex)) + geom_boxplot() + 
    scale_x_discrete(breaks = c("AsstProf", "AssocProf", "Prof"), labels = c("Assistant\nProfessor", 
        "Associate\nProfessor", "Full\nProfessor")) + scale_y_continuous(breaks = c(50000, 
    1e+05, 150000, 2e+05), labels = c("$50K", "$100K", "$150K", "$200K")) + labs(title = "Faculty Salary by Rank and Sex", 
    x = "", y = "")
