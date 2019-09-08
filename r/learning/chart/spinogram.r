library(vcd)

# 棘状图对堆砌条形图进行了重缩放,这样每个条形的高度均为1,每一段的高度即表示比例
Values <- matrix(c(2, 9, 3, 11, 9, 4, 8, 7, 3, 12, 5, 2, 8, 10, 11), nrow = 5, ncol = 3, 
    byrow = T)
spine(Values, main = "Spinogram Example")
# !! spine group by row, bar group by column
