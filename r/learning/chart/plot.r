library(datasets) # Load/unload base packages manually

head(iris)
summary(iris)

?plot # Help for plot()

plot(iris$Species) # Category variable - bar
plot(iris$Petal.Length) # number variable - scatter
plot(iris$Species, iris$Petal.Width) # Category vs number - box
plot(iris$Petal.Length, iris$Petal.Width) # number vs number - scatter
plot(iris) # Entire data frame

# Plot with options
plot(iris$Petal.Length, iris$Petal.Width,
  col = "#cc0000", # Hex code for red
  pch = 19, # Use solid circles for points
  main = "Iris: Petal Length vs. Petal Width",
  xlab = "Petal Length",
  ylab = "Petal Width"
)

# PLOT FORMULAS WITH PLOT() ################################

plot(cos, 0, 2 * pi)
plot(exp, 1, 5)
plot(dnorm, -3, +3)
my_square <- function(x) {
  x * x
}
plot(my_square, -3, 3)

# Formula plot with options
plot(dnorm, -3, +3,
  col = "#cc0000",
  lwd = 5,
  main = "Standard Normal Distribution",
  xlab = "z-scores",
  ylab = "Density"
)

# CLEAN UP #################################################

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear plots
dev.off() # But only if there IS a plot

# Clear console
cat("\014") # ctrl+L