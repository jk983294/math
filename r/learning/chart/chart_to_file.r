v <- c(7, 12, 28, 3, 41)
png(file = "line_chart.png")  # open device
plot(v, type = "o")
dev.off()  # close file

# save to pdf
pdf("my_graph.pdf")
with(mtcars, {
    plot(wt, mpg)
    abline(lm(mpg ~ wt))
    title("Regression of MPG on Weight")
})
dev.off()
