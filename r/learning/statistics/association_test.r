library(vcd)
mytable <- xtabs(~Treatment + Improved, data = Arthritis)
assocstats(mytable)
