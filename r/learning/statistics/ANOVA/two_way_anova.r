# Two way ANOVA
attach(ToothGrowth)
table(supp, dose)
aggregate(len, by = list(supp, dose), FUN = mean)
aggregate(len, by = list(supp, dose), FUN = sd)
dose <- factor(dose)
fit <- aov(len ~ supp * dose)
summary(fit)

# plotting interactions
interaction.plot(dose, supp, len, type = "b", col = c("red", "blue"), pch = c(16, 
    18), main = "Interaction between Dose and Supplement Type")
library(gplots)
plotmeans(len ~ interaction(supp, dose, sep = " "), connect = list(c(1, 3, 5), c(2, 
    4, 6)), col = c("red", "darkgreen"), main = "Interaction Plot with 95% CIs", 
    xlab = "Treatment and Dose Combination")
library(HH)
interaction2wt(len ~ supp * dose)
