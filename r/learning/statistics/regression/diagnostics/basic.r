opar <- par(no.readonly = TRUE)
fit <- lm(weight ~ height, data = women)
par(mfrow = c(2, 2))
plot(fit)  # explanation in readme.md
par(opar)
