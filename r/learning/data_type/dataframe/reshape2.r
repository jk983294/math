library(reshape2)

# input data
mydata <- read.table(header = TRUE, sep = " ", text = "
ID Time X1 X2
1 1 5 6
1 2 3 5
2 1 6 1
2 2 2 4
")

md <- melt(mydata, id = c("ID", "Time"))

# reshaping with aggregation
dcast(md, ID ~ variable, mean)  # group by ID, get mean of X1 and X2
dcast(md, Time ~ variable, mean)  # group by Time, get mean of X1 and X2
dcast(md, ID ~ Time, mean)  # group by (ID, Time), get mean of value

# reshaping without aggregation
dcast(md, ID + Time ~ variable)
dcast(md, ID + variable ~ Time)
dcast(md, ID ~ variable + Time)
