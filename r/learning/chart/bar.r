# bar chart
H <- c(7, 12, 28, 3, 41)
M <- c("Mar", "Apr", "May", "Jun", "Jul")
barplot(H, names.arg = M, xlab = "Month", ylab = "Revenue", col = "blue", main = "Revenue chart", 
    border = "red")

# horizontal bar plot
barplot(H, names.arg = M, xlab = "Month", ylab = "Revenue", col = "blue", main = "Revenue chart", 
    border = "red", horiz = TRUE)

# stacked bar, every bar is accumulated by column in matrix
colors <- c("green", "orange", "brown")
months <- c("Mar", "Apr", "May", "Jun", "Jul")
regions <- c("East", "West", "North")
Values <- matrix(c(2, 9, 3, 11, 9, 4, 8, 7, 3, 12, 5, 2, 8, 10, 11), nrow = 3, ncol = 5, 
    byrow = TRUE)
barplot(Values, main = "total revenue", names.arg = months, xlab = "month", ylab = "revenue", 
    col = colors)
legend("topleft", regions, cex = 1.3, fill = colors)

# group bar chart
barplot(Values, main = "total revenue", names.arg = months, xlab = "month", ylab = "revenue", 
    col = colors, beside = TRUE)
