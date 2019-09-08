library(reshape)
library(MASS)

# Basically, you "melt" data so that each row is a unique id-variable combination. Then you "cast" the melted data into any shape you would like. 

x = data.frame(id = c(1, 1, 2, 2), time = c(1, 2, 1, 2), x1 = c(5, 3, 6, 2), x2 = c(6, 5, 1, 4))

molten = melt(data = x, id=c("id","time"))

# cast(data, formula, function) 
cast(molten, id ~ variable, sum) 
cast(molten, id ~ variable, mean)
cast(molten, time ~ variable, min)
cast(molten, time ~ variable, max)
cast(molten, id + time ~ variable, sum) 
cast(molten, id + time ~ variable, mean)
