x <- c(-11.2, -11.3, 9.4, 13.2, 21, 8.1, 13.6, 2.2, 12.2, 4.1, 13.1, 3.3, 19)
hist(x)
hist(x, xlim = c(-15, 25), breaks = 4, col = "red")

# 轴须图
hist(x,
    freq = FALSE, # Axis shows density, not freq.
    col = "thistle1", # Color for histogram
    main = "Histogram of x",
    xlab = "x"
)

# Add a normal distribution
curve(dnorm(x, mean = mean(x), sd = sd(x)),
    col = "thistle4", # Color of curve
    lwd = 2, # Line width of 2 pixels
    add = TRUE
) # Superimpose on previous graph

# Add kernel density estimators
lines(density(x), col = "blue", lwd = 2)

# Add a rug plot
rug(x, lwd = 2, col = "gray")