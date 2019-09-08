attach(mtcars)
r <- sqrt(disp/pi)
symbols(wt, mpg, circle = r, inches = 0.3, fg = "white", bg = "lightblue", main = "Bubble Plot with point size proportional to displacement", 
    ylab = "Miles Per Gallon", xlab = "Weight of Car (lbs/1000)")
text(wt, mpg, rownames(mtcars), cex = 0.6)
detach(mtcars)
