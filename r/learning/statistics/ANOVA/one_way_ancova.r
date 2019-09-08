data(litter, package = "multcomp")
attach(litter)
table(dose)
aggregate(weight, by = list(dose), FUN = mean)
fit <- aov(weight ~ gesttime + dose)
summary(fit)

# Visualizing a one-way ANCOVA
library(HH)
ancova(weight ~ gesttime + dose, data = litter)
