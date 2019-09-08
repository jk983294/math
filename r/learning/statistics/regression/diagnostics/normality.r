# assess normality
library(car)
states <- as.data.frame(state.x77[, c("Murder", "Population", "Illiteracy", "Income", 
    "Frost")])
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)
qqPlot(fit, labels = row.names(states), id.method = "identify", simulate = TRUE, 
    main = "Q-Q Plot")
