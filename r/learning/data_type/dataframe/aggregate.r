options(digits = 3)
attach(mtcars)
# group by (cyl, gear) and Calculate mean of each group for left columns
aggregate(mtcars, by = list(Group.cyl = cyl, Group.gear = gear), FUN = mean, na.rm = TRUE)
