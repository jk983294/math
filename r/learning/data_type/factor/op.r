# facor stores the vector along with the distinct values of the elements in the
# vector as labels
apple_colors <- c("green", "green", "yellow", "red", "red", "red", "green")
factor_apple <- factor(apple_colors)
nlevels(factor_apple)  # 3
levels(factor_apple) # "green"  "red"    "yellow"
