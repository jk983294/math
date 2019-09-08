library(vcd)
# one way table
mytable <- table(Arthritis$Improved)  # 频数统计表 on column Improved
prop.table(mytable)  # proportions
prop.table(mytable) * 100  # percentages

# two way table
mytable <- xtabs(~Treatment + Improved, data = Arthritis)
mytable  # frequencies
# 边际频数和比例
margin.table(mytable, 1)  # row sums
margin.table(mytable, 2)  # column sums
prop.table(mytable)  # cell proportions
prop.table(mytable, 1)  # row proportions
prop.table(mytable, 2)  # column proportions
addmargins(mytable)  # add row and column sums to table

# more complex tables
addmargins(prop.table(mytable))
addmargins(prop.table(mytable, 1), 2)
addmargins(prop.table(mytable, 2), 1)

# Two way table using CrossTable
library(gmodels)
CrossTable(Arthritis$Treatment, Arthritis$Improved)

# high dimension table
mytable <- xtabs(~Treatment + Sex + Improved, data = Arthritis)
ftable(mytable)
ftable(prop.table(mytable, c(1, 2)))
ftable(addmargins(prop.table(mytable, c(1, 2)), 3))
